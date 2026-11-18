# Detecção e classificação de câncer por metilação do DNA

Trabalho da disciplina **Genômica Computacional**. Baseado em
Liu, M. C. et al. (2020), *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*, Annals of Oncology 31(6):745–759.

O artigo usa dados privados (estudos CCGA/STRIVE, da GRAIL). Aqui construímos um classificador **Random Forest** com dados públicos:

- **Experimento principal — tecido tumoral (TCGA, 450K):** 6 tipos de câncer (cólon, estômago, pulmão, fígado, mama, rim), 1.290 amostras. Etapa 1: tumor vs normal; etapa 2: tecido de origem. Validação por *nested CV* agrupada por paciente, com auditoria de vazamento.
- **Tentativa documentada — plasma (GEO GSE186575, cfMeDIP-seq):** colorretal, gástrico e saudáveis. A detecção ficou no nível do acaso, dominada por lote e idade; serve de contraste com o artigo.

## Estrutura

```
├── PLANO_EXECUCAO.md        # plano, status das etapas e log de decisões
├── docs/                    # enunciado e artigo (PDFs não versionados: material da disciplina/editora)
├── notebooks/
│   ├── trabalho_final.ipynb # notebook de ENTREGA (roda no Colab)
│   └── exploracao/          # notebooks de teste/rascunho
├── data/
│   ├── raw/                 # GEO (versionado, ~22 MB); raw/tcga/ baixado pelo notebook (~2,7 GB, fora do git)
│   ├── reference/           # anotação gene → cromossomo (Ensembl hg38)
│   └── processed/           # matrizes intermediárias geradas pelo notebook
├── results/
│   ├── figures/             # gráficos gerados
│   └── tables/              # métricas e tabelas geradas
├── relatorio/               # relatório técnico
└── scripts/                 # scripts auxiliares (ex.: gerar notebooks)
```

## Como rodar

**Google Colab:** abrir `notebooks/trabalho_final.ipynb` direto do GitHub (Arquivo → Abrir notebook → GitHub) e executar tudo. Não precisa subir nada:
- por padrão, os dados já tratados (~860 MB) são baixados de uma pasta pública do Google Drive via `gdown`, sem login;
- se o Drive falhar (ou com `USAR_DRIVE = False`), tudo é baixado das fontes originais (TCGA via UCSC Xena, ~2,7 GB; GEO; Ensembl) e processado do zero.

Tempo estimado: ~30–50 min no Colab gratuito (pico de memória ~6 GB). A auditoria pode ser desligada com `RODAR_AUDITORIA = False`.

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
| TCGA 450K: COAD, STAD, LUAD, LIHC, BRCA, KIRC + probe map | [UCSC Xena](https://xenabrowser.net/datapages/) | Público |
| GSE186575 (cfMeDIP-seq de plasma: CRC, gástrico, saudáveis) | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186575) | Público |
| Anotação de genes (hg38) | [Ensembl BioMart](https://www.ensembl.org/biomart) | Público |
| Dados já tratados (cópia para agilizar o Colab) | [Google Drive](https://drive.google.com/drive/folders/1PKbiUFQM2cswK3ShJszUrZzyWwDokuBR) | Público (leitura) |

Artigo de referência: Liu et al. (2020), doi:[10.1016/j.annonc.2020.02.011](https://doi.org/10.1016/j.annonc.2020.02.011).
