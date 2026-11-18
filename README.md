# Detecção e classificação de câncer por metilação do DNA

Trabalho da disciplina **Genômica Computacional**. Baseado em
Liu, M. C. et al. (2020), *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*, Annals of Oncology 31(6):745–759.

O artigo usa dados privados (estudos CCGA/STRIVE, da GRAIL). Aqui construímos um classificador **Random Forest** com dados públicos de cfDNA plasmático (**GEO GSE186575**): câncer colorretal, gástrico e controles saudáveis.

## Estrutura

```
├── PLANO_EXECUCAO.md        # plano, status das etapas e log de decisões
├── docs/                    # enunciado e artigo de referência (PDF)
├── notebooks/
│   ├── trabalho_final.ipynb # notebook de ENTREGA (roda no Colab)
│   └── exploracao/          # notebooks de teste/rascunho
├── data/
│   ├── raw/                 # dados originais do GEO (versionados, ~22 MB)
│   ├── reference/           # anotação gene → cromossomo (Ensembl hg38)
│   └── processed/           # matrizes intermediárias geradas pelo notebook
├── results/
│   ├── figures/             # gráficos gerados
│   └── tables/              # métricas e tabelas geradas
├── relatorio/               # relatório técnico
└── scripts/                 # scripts auxiliares (ex.: gerar notebooks)
```

## Como rodar

**Google Colab:** abrir `notebooks/trabalho_final.ipynb` direto do GitHub (Arquivo → Abrir notebook → GitHub). A primeira célula clona o repositório, então os dados já vêm junto. Se faltar algum arquivo, ele é baixado do FTP público do GEO.

**Local:**
```bash
python -m venv .venv
.venv/Scripts/activate        # Windows  (Linux/Mac: source .venv/bin/activate)
pip install -r requirements.txt
jupyter notebook notebooks/trabalho_final.ipynb
```

## Fontes de dados

| Dado | Origem | Acesso |
|---|---|---|
| GSE186575 (cfMeDIP-seq de plasma: CRC, gástrico, saudáveis) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186575) | Público |
| Anotação de genes (hg38) | [Ensembl BioMart](https://www.ensembl.org/biomart) | Público |
