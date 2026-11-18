"""Gera o notebook de modelagem (rascunho) com Random Forest no GSE186575.

Uso: python scripts/build_model_notebook.py
"""
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "exploracao" / "02_modelagem_rf.ipynb"

md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell

cells = [
    md("""# 02 — Modelagem com Random Forest (notebook de testes)

Duas tarefas, espelhando o artigo (Liu et al. 2020):
1. **Detecção:** câncer (CRC + STAD) vs saudável (NC)
2. **Tecido de origem (TOO):** CRC vs STAD, só nas amostras com câncer
3. **Fluxo hierárquico:** aplica (1) e, nos positivos, (2)

Regras anti-vazamento:
- Toda etapa que **olha os dados de várias amostras** (filtro de promotores, seleção de features) fica **dentro do `Pipeline`** → é ajustada só no treino de cada fold.
- `log2(CPM+1)` é aplicado amostra a amostra (não usa estatística de outras amostras) → pode ficar fora.
- Cada amostra do GEO é um indivíduo diferente (um GSM por pessoa) → não há risco de o mesmo paciente cair em treino e teste."""),
    code("""import sys, gzip, urllib.request
from pathlib import Path

IN_COLAB = "google.colab" in sys.modules
REPO_URL = "https://github.com/<usuario>/GenomicaComputacional.git"  # TODO: ajustar após criar o repo
if IN_COLAB and not Path("GenomicaComputacional").exists():
    !git clone -q {REPO_URL}

def find_root():
    for p in [Path.cwd(), *Path.cwd().parents, Path("GenomicaComputacional").resolve()]:
        if (p / "PLANO_EXECUCAO.md").exists():
            return p
    return Path.cwd()

ROOT = find_root()
RAW, REF = ROOT / "data" / "raw", ROOT / "data" / "reference"
FIG, TAB = ROOT / "results" / "figures", ROOT / "results" / "tables"
for d in (RAW, FIG, TAB): d.mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import (StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV,
                                     train_test_split, cross_val_predict)
from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay, classification_report,
                             roc_auc_score, roc_curve, f1_score, precision_score, recall_score, accuracy_score)
from sklearn.inspection import permutation_importance

sns.set_theme(style="whitegrid")
RANDOM_STATE = 42"""),
    md("## 1. Dados (mesmo pré-processamento validado na EDA)"),
    code("""GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186575"
FILES = {"GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz": f"{GEO_FTP}/suppl/GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz",
         "GSE186575_series_matrix.txt.gz": f"{GEO_FTP}/matrix/GSE186575_series_matrix.txt.gz"}
for name, url in FILES.items():
    if not (RAW / name).exists():
        urllib.request.urlretrieve(url, RAW / name)

def read_geo_meta(path):
    rows = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\\n").split("\\t")
                rows.setdefault(k, []).append([x.strip('"') for x in v])
    meta = pd.DataFrame({"title": rows["!Sample_title"][0], "gsm": rows["!Sample_geo_accession"][0]})
    for ch in rows["!Sample_characteristics_ch1"]:
        meta[ch[0].split(":")[0].strip().lower()] = [x.split(":", 1)[1].strip() for x in ch]
    return meta

meta = read_geo_meta(RAW / "GSE186575_series_matrix.txt.gz")
meta["grupo"] = meta["title"].str.split("-").str[0]
meta["lote_mix"] = meta["title"].str.contains("mix")
meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
meta = meta.set_index("title")

X_cpm = pd.read_csv(RAW / "GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz", sep="\\t", index_col=0).T.loc[meta.index]

# anotação e remoção de chrX/chrY/MT/sem anotação (decisão D6)
ann = pd.read_csv(REF / "ensembl_genes_hg38.tsv", sep="\\t").drop_duplicates("Gene stable ID").set_index("Gene stable ID")
gene_id = X_cpm.columns.str.replace("promoter_", "").str.split("|").str[0].str.split(".").str[0]
feat = pd.DataFrame({"gene_id": gene_id,
                     "gene": ann["Gene name"].reindex(gene_id).values,
                     "chrom": ann["Chromosome/scaffold name"].reindex(gene_id).values}, index=X_cpm.columns)
auto = feat.index[feat["chrom"].isin([str(i) for i in range(1, 23)])]

X = np.log2(X_cpm[auto] + 1)   # transformação amostra a amostra → sem vazamento
print("X:", X.shape)
meta["grupo"].value_counts()"""),
    md("""## 2. Pipeline

1. **`PresenceFilter`** — remove promotores com sinal (CPM > 0) em menos de 50% das amostras **de treino** (não supervisionado, mas ajustado só no treino).
2. **`SelectKBest(f_classif)`** — teste F (ANOVA) entre classes; mantém os *k* promotores mais discriminativos (**supervisionado** → obrigatoriamente dentro da validação).
3. **`RandomForestClassifier`** — `class_weight="balanced"` compensa o desbalanceamento.

Hiperparâmetros (*k* e `max_features`) escolhidos por **GridSearchCV interno**; desempenho estimado por **CV externa** (validação aninhada / *nested CV*)."""),
    code("""class PresenceFilter(BaseEstimator, TransformerMixin):
    \"\"\"Mantém features com valor > 0 em pelo menos `min_frac` das amostras de treino.\"\"\"
    def __init__(self, min_frac=0.5):
        self.min_frac = min_frac
    def fit(self, X, y=None):
        self.mask_ = (np.asarray(X) > 0).mean(axis=0) >= self.min_frac
        return self
    def transform(self, X):
        return np.asarray(X)[:, self.mask_]
    def get_support(self):
        return self.mask_

def make_pipeline():
    return Pipeline([
        ("filtro", PresenceFilter(min_frac=0.5)),
        ("selecao", SelectKBest(f_classif)),
        ("rf", RandomForestClassifier(n_estimators=500, class_weight="balanced",
                                      random_state=RANDOM_STATE, n_jobs=-1)),
    ])

PARAM_GRID = {"selecao__k": [100, 500, 2000], "rf__max_features": ["sqrt", 0.1]}

def make_search():
    inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    return GridSearchCV(make_pipeline(), PARAM_GRID, cv=inner, scoring="roc_auc", n_jobs=1)

def feature_names(pipe, cols):
    \"\"\"Recupera os nomes dos promotores que sobreviveram ao filtro + seleção.\"\"\"
    cols = np.asarray(cols)[pipe.named_steps["filtro"].get_support()]
    return cols[pipe.named_steps["selecao"].get_support()]"""),
    md("""## 3. Funções de avaliação

- **Hold-out (70/30 estratificado):** um conjunto de teste fixo, nunca visto no treino nem na escolha de hiperparâmetros → matriz de confusão e curva ROC "apresentáveis".
- **CV externa repetida (5 folds × 3 repetições):** com só 98 amostras, um único split é instável; a CV repetida dá média ± desvio."""),
    code("""def metricas(y_true, y_pred, y_prob, pos_label=1):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {"acuracia": accuracy_score(y_true, y_pred),
            "precisao": precision_score(y_true, y_pred, zero_division=0),
            "sensibilidade": recall_score(y_true, y_pred),
            "especificidade": tn / (tn + fp),
            "f1": f1_score(y_true, y_pred),
            "roc_auc": roc_auc_score(y_true, y_prob)}

def nested_cv(Xt, yt, n_repeats=3):
    outer = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=RANDOM_STATE)
    res, params = [], []
    for tr, te in outer.split(Xt, yt):
        s = make_search().fit(Xt.iloc[tr], yt.iloc[tr])
        prob = s.predict_proba(Xt.iloc[te])[:, 1]
        res.append(metricas(yt.iloc[te], (prob >= .5).astype(int), prob))
        params.append(s.best_params_)
    res = pd.DataFrame(res)
    print("Nested CV (média ± desvio):")
    display(res.agg(["mean", "std"]).T.round(3))
    print("Hiperparâmetros escolhidos nos folds:")
    display(pd.DataFrame(params).value_counts().rename("vezes"))
    return res

def holdout(Xt, yt, nomes, tag):
    X_tr, X_te, y_tr, y_te = train_test_split(Xt, yt, test_size=0.3, stratify=yt, random_state=RANDOM_STATE)
    s = make_search().fit(X_tr, y_tr)
    prob = s.predict_proba(X_te)[:, 1]; pred = (prob >= .5).astype(int)
    print("melhores hiperparâmetros:", s.best_params_)
    m = metricas(y_te, pred, prob); display(pd.Series(m).round(3).to_frame(f"teste ({len(y_te)} amostras)"))
    print(classification_report(y_te, pred, target_names=nomes, zero_division=0))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ConfusionMatrixDisplay(confusion_matrix(y_te, pred), display_labels=nomes).plot(ax=ax[0], cmap="Blues", colorbar=False)
    ax[0].set_title(f"{tag}: matriz de confusão (teste)")
    fpr, tpr, _ = roc_curve(y_te, prob)
    ax[1].plot(fpr, tpr, lw=2, label=f"AUC = {m['roc_auc']:.3f}"); ax[1].plot([0, 1], [0, 1], "k--", lw=1)
    ax[1].set(xlabel="1 − especificidade", ylabel="sensibilidade", title=f"{tag}: curva ROC (teste)"); ax[1].legend()
    plt.tight_layout(); plt.savefig(FIG / f"{tag}_holdout.png", dpi=150); plt.show()
    return s.best_estimator_, X_tr, X_te, y_tr, y_te"""),
    md("""## 4. Tarefa 1 — Detecção: câncer vs saudável"""),
    code("""y_det = (meta["grupo"] != "NC").astype(int)   # 1 = câncer, 0 = saudável
print(y_det.value_counts().rename({1: "câncer", 0: "saudável"}))
res_det = nested_cv(X, y_det)"""),
    code("""modelo_det, Xd_tr, Xd_te, yd_tr, yd_te = holdout(X, y_det, ["saudável", "câncer"], "deteccao")"""),
    md("""## 5. Tarefa 2 — Tecido de origem: CRC vs STAD (só amostras com câncer)"""),
    code("""cancer = meta.index[meta["grupo"] != "NC"]
X_too = X.loc[cancer]
y_too = (meta.loc[cancer, "grupo"] == "STAD").astype(int)   # 1 = STAD (estômago), 0 = CRC (colorretal)
print(y_too.value_counts().rename({1: "STAD", 0: "CRC"}))
res_too = nested_cv(X_too, y_too)"""),
    code("""modelo_too, Xt_tr, Xt_te, yt_tr, yt_te = holdout(X_too, y_too, ["CRC", "STAD"], "too")"""),
    md("""## 6. Fluxo hierárquico (detecção → TOO), como no artigo

Usa as predições **fora da amostra** (*out-of-fold*, via `cross_val_predict`) para todas as 98 amostras:
- etapa 1 decide câncer vs saudável;
- etapa 2 (treinada só com câncer) atribui CRC/STAD a quem foi positivo.

Hiperparâmetros fixados nos valores mais escolhidos na nested CV para manter o custo baixo."""),
    code("""def pipe_fixo(res_search_params):
    p = make_pipeline(); p.set_params(**res_search_params); return p

cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
p_det = pipe_fixo({"selecao__k": 500, "rf__max_features": "sqrt"})
p_too = pipe_fixo({"selecao__k": 500, "rf__max_features": "sqrt"})

oof_det = pd.Series(cross_val_predict(p_det, X, y_det, cv=cv5, method="predict_proba")[:, 1], index=X.index)
oof_too = pd.Series(cross_val_predict(p_too, X_too, y_too, cv=cv5, method="predict_proba")[:, 1], index=X_too.index)

# Para saudáveis classificados como câncer, a etapa 2 precisa de um modelo treinado em todos os cânceres
p_too_full = pipe_fixo({"selecao__k": 500, "rf__max_features": "sqrt"}).fit(X_too, y_too)
final = pd.Series("NC", index=X.index)
pos = oof_det.index[oof_det >= .5]
for s in pos:
    p = oof_too[s] if s in oof_too.index else p_too_full.predict_proba(X.loc[[s]])[0, 1]
    final[s] = "STAD" if p >= .5 else "CRC"

ordem = ["NC", "CRC", "STAD"]
cm = confusion_matrix(meta["grupo"], final, labels=ordem)
ConfusionMatrixDisplay(cm, display_labels=ordem).plot(cmap="Blues", colorbar=False)
plt.title("Fluxo hierárquico (out-of-fold, 98 amostras)"); plt.savefig(FIG / "hierarquico_cm.png", dpi=150, bbox_inches="tight"); plt.show()
print(classification_report(meta["grupo"], final, labels=ordem, zero_division=0))"""),
    md("""## 7. Confundidores: lote `mix` e idade

**7a. Amostras que cruzam o lote.** Se o modelo de detecção aprendeu "lote mix ⇒ saudável", ele deve errar justamente:
- nos **6 saudáveis não-mix** (chamando-os de câncer), e
- nos **2 CRC mix** (chamando-os de saudáveis)."""),
    code("""tab = meta[["grupo", "lote_mix", "age"]].copy()
tab["prob_cancer"] = oof_det.round(3)
tab["pred"] = np.where(oof_det >= .5, "câncer", "saudável")
tab["acertou"] = (oof_det >= .5).astype(int) == y_det
display(tab.groupby(["grupo", "lote_mix"]).agg(n=("acertou", "size"), taxa_acerto=("acertou", "mean"),
                                               prob_cancer_media=("prob_cancer", "mean")).round(3))
display(tab[(tab.grupo == "NC") & ~tab.lote_mix | (tab.grupo == "CRC") & tab.lote_mix])"""),
    md("""**7b. A idade sozinha prevê câncer?** Um RF usando **só a idade** serve de linha de base: se ele já fosse muito bom, parte do desempenho do modelo de metilação poderia vir da idade."""),
    code("""rf_idade = RandomForestClassifier(n_estimators=500, min_samples_leaf=5, class_weight="balanced", random_state=RANDOM_STATE)
prob_idade = cross_val_predict(rf_idade, meta[["age"]], y_det, cv=cv5, method="predict_proba")[:, 1]
print(f"AUC só com idade: {roc_auc_score(y_det, prob_idade):.3f}")
print(f"AUC metilação (out-of-fold): {roc_auc_score(y_det, oof_det):.3f}")"""),
    md("""## 8. Quais promotores mais contribuem?

*Permutation importance* no conjunto de **teste** do hold-out: embaralha uma feature por vez e mede quanto o AUC cai. É mais confiável que a importância por impureza (que favorece features com muitos valores distintos)."""),
    code("""def importancias(modelo, X_te, y_te, cols, tag, top=20):
    nomes = feature_names(modelo, cols)
    Xs = modelo[:-1].transform(X_te)                       # aplica filtro + seleção
    pi = permutation_importance(modelo.named_steps["rf"], Xs, y_te, scoring="roc_auc",
                                n_repeats=20, random_state=RANDOM_STATE, n_jobs=-1)
    imp = (pd.DataFrame({"promotor": nomes, "imp_permutacao": pi.importances_mean,
                         "imp_impureza": modelo.named_steps["rf"].feature_importances_})
             .join(feat, on="promotor").sort_values("imp_permutacao", ascending=False))
    imp.to_csv(TAB / f"importancias_{tag}.csv", index=False)
    t = imp.head(top)
    plt.figure(figsize=(7, 6)); sns.barplot(data=t, y="gene", x="imp_permutacao", color="C0")
    plt.title(f"{tag}: top {top} promotores (permutation importance)"); plt.tight_layout()
    plt.savefig(FIG / f"importancias_{tag}.png", dpi=150); plt.show()
    return imp

imp_det = importancias(modelo_det, Xd_te, yd_te, X.columns, "deteccao")
imp_det.head(15)"""),
    code("""imp_too = importancias(modelo_too, Xt_te, yt_te, X_too.columns, "too")
imp_too.head(15)"""),
    md("""## 9. Demonstração de vazamento (para o relatório)

Mesmo modelo, mas com a seleção de features feita **no dataset inteiro antes** da validação cruzada (errado) vs **dentro** de cada fold (certo).
Para deixar o efeito evidente, usamos **rótulos embaralhados** (aleatórios): não existe sinal real, então o AUC honesto deveria ficar perto de 0,5."""),
    code("""rng = np.random.default_rng(RANDOM_STATE)
y_rand = pd.Series(rng.permutation(y_det.values), index=y_det.index)

# ERRADO: seleciona os 500 promotores mais associados ao rótulo usando TODAS as amostras
sel = SelectKBest(f_classif, k=500).fit(X, y_rand)
rf = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1)
prob_errado = cross_val_predict(rf, X.loc[:, sel.get_support()], y_rand, cv=cv5, method="predict_proba")[:, 1]

# CERTO: seleção dentro do pipeline
prob_certo = cross_val_predict(pipe_fixo({"selecao__k": 500, "rf__max_features": "sqrt"}),
                               X, y_rand, cv=cv5, method="predict_proba")[:, 1]
print(f"Rótulos aleatórios — AUC com vazamento: {roc_auc_score(y_rand, prob_errado):.3f}")
print(f"Rótulos aleatórios — AUC sem vazamento: {roc_auc_score(y_rand, prob_certo):.3f}")"""),
    code("""# salva um resumo das métricas
resumo = pd.concat({"deteccao": res_det.agg(["mean", "std"]).T, "too": res_too.agg(["mean", "std"]).T}, axis=0).round(3)
resumo.to_csv(TAB / "metricas_nested_cv.csv"); resumo"""),
]

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
    "language_info": {"name": "python"},
})
nbf.write(nb, OUT)
print("notebook escrito em", OUT)
