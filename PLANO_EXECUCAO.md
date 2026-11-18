# Plano de Execução — Detecção/classificação de câncer por metilação do DNA

**Disciplina:** Genômica Computacional
**Artigo de referência:** Liu, M. C. et al. (2020). *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*. Annals of Oncology, 31(6), 745–759. DOI: 10.1016/j.annonc.2020.02.011
**Documentos-fonte:** [docs/atividade_genomica_methylacao_multicancer.pdf](docs/atividade_genomica_methylacao_multicancer.pdf), [docs/artigo_referencia.pdf](docs/artigo_referencia.pdf)

> Legenda de status: `[ ]` pendente · `[~]` em andamento · `[x]` concluído
> **Próxima etapa:** ver a seção [Status atual](#status-atual) no fim do arquivo.

---

## 0. Entregáveis (o que precisa existir no final)

- [ ] **Relatório técnico** — síntese do artigo, investigação dos datasets, justificativa da escolha, metodologia, resultados, discussão, referências
- [ ] **Notebook Google Colab reproduzível** — leitura dos dados → preparação → treino → avaliação
- [ ] **Registro explícito** de: fontes de dados, definição das classes, estratégia de validação (pode ficar no relatório + início do notebook)

Critério central do professor: **escolha de dados adequada** e **prevenção de vazamento de informação** (data leakage).

---

## 1. Parte I — Estudo do artigo

### 1.1 Fundamentação biológica
- [ ] Metilação do DNA (5mC em dinucleotídeos CpG, ilhas CpG, promotores) e sua relação com câncer (hipermetilação de promotores de genes supressores, hipometilação global)
- [ ] cfDNA: fragmentos de DNA liberados no plasma (apoptose/necrose); ctDNA = fração tumoral; uso em biópsia líquida
- [ ] Por que metilação detecta câncer **e** tecido de origem (TOO): padrão de metilação é específico de tipo celular/tecido → assinatura do tecido que liberou o DNA

### 1.2 Caracterização das amostras (já extraído do artigo)
| Conjunto | Câncer | Não-câncer | Total |
|---|---|---|---|
| Treinamento | 1.742 | 2.978 (1.391 CCGA + 1.587 STRIVE) | 4.720 |
| Validação | 740 | 1.229 (614 CCGA + 615 STRIVE) | 1.969 |
| Total | 2.482 | 4.207 | 6.689 |

- Analisáveis na análise primária: 4.316 (treino 3.052; validação 1.264)
- 354 amostras do CCGA reservadas para um *tumor biopsy reference set*

- [ ] **CCGA (NCT02889978):** estudo prospectivo, multicêntrico, caso-controle, observacional, 15.254 participantes de 142 centros na América do Norte. Dividido em 3 sub-estudos: (1) descoberta — comparou WGBS, WGS e sequenciamento direcionado; (2) treino/validação do ensaio de metilação direcionada (**este artigo**); (3) validação adicional. Fornece os casos de câncer e parte dos controles.
- [ ] **STRIVE (NCT03085888):** coorte prospectiva de mulheres em mamografia de rastreio (99.286 participantes, 35 centros). Forneceu **controles não-câncer extras** (de um único centro) para treinar e medir especificidade >99% com confiança >90%.
- [ ] **Tecnologia:** cfDNA plasmático (até 75 ng) → conversão por bissulfito (Zymo EZ-96) → biblioteca dual-index → captura híbrida com painel customizado (Twist) → sequenciamento 150 bp paired-end, Illumina NovaSeq. Profundidade mediana 113 M fragmentos (139X único on-target).
- [ ] **Escolha das regiões:** baseada em WGBS do 1º sub-estudo (cfDNA, biópsias FFPE, leucócitos, tecidos comerciais) + arrays de metilação do TCGA → regiões do hg19 com padrões câncer/tecido-específicos em cfDNA vs controles → algoritmo customizado → 103.456 regiões (17,2 Mb), ~1,1 M CpGs.
- [ ] **Separação dos dados e vazamento:** treino e validação independentes; classificador **travado** antes de liberar a validação; análises **duplo-cegas**; equipe de integridade de dados separada; distribuição de tipos de câncer/estádio pré-especificada e balanceada entre centros; análise *post hoc* de balanceamento por centro.
  - Pontos críticos para discutir: (a) regiões do painel foram escolhidas com dados de descoberta (1º sub-estudo) — se houvesse sobreposição com validação seria vazamento; (b) controles STRIVE são só mulheres e de um único centro → possível viés de lote/sexo/centro (artigo mostra FPR similar CCGA vs STRIVE: 0,7% vs 0,6%); (c) modelos treinados em 11.154 amostras de 5.854 participantes **com CV**, garantindo que folds held-out não fossem vistos — note: múltiplas amostras por participante → separação por indivíduo é necessária.

### 1.3 Estratégia de classificação
- [ ] Problema (i): **detecção** câncer vs não-câncer (binário), operando em especificidade fixa alta (~99,3%)
- [ ] Problema (ii): **TOO** — só nas amostras positivas, classificar o tecido de origem (multiclasse)
- [ ] Arquitetura: *source models* por região (reconhecem padrões de metilação semelhantes aos de cada tipo de câncer) → features → **par de regressões logísticas em ensemble** (uma para câncer/não-câncer, outra para TOO). Thresholds também ajustados com CV.
- [ ] Resultados-chave: especificidade 99,3%; sensibilidade estádio I–III 67,3% (12 cânceres pré-especificados) e 43,9% (todos); sensibilidade cresce com estádio (18% → 93%); TOO previsto em 96% dos positivos, acurácia 93%.

---

## 2. Parte II — Busca e avaliação de datasets públicos

### 2.1 CCGA / STRIVE — dados não disponíveis publicamente
- [x] Artigo **não tem** seção de *data availability* nem accession (EGA/dbGaP/GEO) — dados são da GRAIL (financiadora).
- [x] ClinicalTrials.gov API: **STRIVE (NCT03085888) → `ipdSharing: NO`** (sem compartilhamento de dados individuais). **CCGA (NCT02889978)** → nenhum plano de compartilhamento declarado no registro.
- Conclusão: **não é possível obter** os dados de metilação do CCGA/STRIVE → justificar uso de dataset alternativo.
- Exemplo de **acesso controlado** encontrado: compêndio pan-câncer de cfMeDIP-seq (Zeng et al., *Nature Cancer* 2026, 1.294 amostras, 11 tipos de câncer) — dados brutos só via EGA com *Data Access Agreement* (UHN DAC) ou pedido aos autores. Inviável no prazo; bom exemplo para discutir “público vs controlado”.

### 2.2 Tabela comparativa dos candidatos (verificada em 2026-11-17 via API do GEO/Xena)

| Campo | **GSE186575** | **GSE81314** | **TCGA (UCSC Xena)** | GSE314062 | GSE275598 | GSE122126 |
|---|---|---|---|---|---|---|
| Link | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186575) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81314) | [Xena](https://xenabrowser.net/datapages/?hub=https://tcga.xenahubs.net:443) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314062) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275598) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122126) |
| Publicação | PMID 37992683 (multi-ômica plasmática, câncer GI, China) | Li et al. 2017, *Cell Research* (5hmC em cfDNA) | TCGA / PanCanAtlas | PMID 42524278 | PMID 39187719 | Moss et al. 2018, *Nat Commun* (PMID 30498206) |
| Nº amostras | 98 | 126 GSMs → **matriz com 68 colunas** (64 pacientes úteis) | ~400–900 por coorte (ver 2.3) | 5 | 93 | 101 (59 cfDNA) |
| Classes | CRC 32 · gástrico (STAD) 35 · saudável 31 | Pulmão 15 · HCC 10 · pâncreas 7 · gástrico 5 · cólon 4 · mama 4 · GBM 4 · HBV 7 · saudável 8 | Tumor primário vs normal adjacente, por tipo de câncer | 3 pacientes HNC, sem controles | 7 subtipos de tumor cerebral pediátrico + 16 não-neoplásicos | Normal, sepse, transplante de ilhota, poucos cânceres (4 cólon, 4 pulmão, 4 CUP) |
| Controles | 31 doadores saudáveis | 8 saudáveis + 7 HBV (não-câncer) | Tecido normal adjacente (poucos em alguns tipos) | Nenhum | 16 não-neoplásicos | Saudáveis |
| Material | **Plasma (cfDNA)** | **Plasma (cfDNA)** | **Tecido tumoral** | Plasma | Plasma | Plasma + células isoladas |
| Tecnologia | cfMeDIP-seq (imunoprecipitação de DNA metilado) | hMe-Seal (5-**hidroxi**metilcitosina, não 5mC) | Illumina 450K array (beta values) | Bisulfito direcionado (11 marcadores) | cfMeDIP-seq | 450K/EPIC array |
| Formato | Matriz CPM por **promotor** (58.288 promotores × 98), 22 MB | Matriz por **gene** (20.244 genes × 68), 9 MB | Matriz beta CpG × amostra (~485k probes), 100–800 MB/coorte | RAW tar | CSV CPM normalizado | RAW tar (IDATs) |
| Acesso | Público | Público | Público | Público | Público | Público |
| Viabilidade | ✅ **Alta** — plasma + controles + 2 tipos (permite detecção **e** TOO) | ⚠️ Média — multiclasse, mas n muito pequeno por classe e sinal é 5hmC | ✅ Alta para **TOO em tecido**; não é cfDNA | ❌ n=5, sem controles | ⚠️ Pediátrico/cérebro, foge do tema multicâncer adulto | ❌ Poucos cânceres; feito para deconvolução |

**Observações importantes encontradas na inspeção:**
- **GSE186575:** 25 dos 31 controles têm título `NC-PKU-mix*` e só 2 CRC têm `mix` → provável **lote (batch) diferente** entre casos e controles. Risco de o classificador aprender o lote e não o câncer (*confounding*). Registrar e investigar no EDA (PCA colorida por lote).
- **GSE81314:** a matriz inclui colunas que **não são pacientes** (`Blood`, `Blood-input`, `Input` = DNA genômico e controles técnicos) → remover. Alguns pacientes têm 2 GSMs (réplicas/RNA), mas a matriz já está por paciente.
- **TCGA:** o arquivo PanCan único tem **14 GB** (inviável no Colab) → usar arquivos por coorte. STAD tem só 2 normais.

### 2.3 TCGA — amostras por coorte (450K, Xena)
| Coorte | Tumor | Normal | Tamanho |
|---|---|---|---|
| BRCA (mama) | 790 | 98 | 783 MB |
| LUAD (pulmão adeno) | 460 | 32 | 450 MB |
| LUSC (pulmão escamoso) | 372 | 43 | 372 MB |
| COAD (cólon) | 299 | 38 | 309 MB |
| READ (reto) | 99 | 7 | 99 MB |
| LIHC (fígado) | 379 | 50 | 386 MB |
| KIRC (rim) | 320 | 160 | 423 MB |
| STAD (estômago) | 396 | 2 | 360 MB |
| PAAD (pâncreas) | 185 | 10 | 178 MB |
| PRAD (próstata) | 499 | 50 | 483 MB |
| HNSC (cabeça/pescoço) | 530 | 50 | 516 MB |
| THCA (tireoide) | 515 | 56 | 493 MB |
| UCEC (endométrio) | 432 | 46 | 429 MB |
| ESCA (esôfago) | 186 | 16 | 186 MB |

> Lembrete do enunciado: classificar **tecido tumoral (TCGA)** ≠ detectar câncer em **cfDNA plasmático**. A definição do problema deve seguir o tipo de dado escolhido.

### 2.4 Proposta de uso (aguardando aprovação)
- **Experimento principal — GSE186575 (plasma):** fluxo hierárquico igual ao do artigo, em escala pequena:
  1. Detecção: câncer (CRC + STAD) vs saudável → binário
  2. TOO: CRC vs STAD, só nas amostras com câncer
- **Experimento complementar (opcional) — TCGA tecido:** multiclasse de tecido de origem com mais amostras (ex.: COAD, STAD, LUAD, LIHC, PAAD, BRCA), deixando claro que é **tecido**, não cfDNA.
- GSE81314 fica como alternativa/menção na discussão (5hmC ≠ 5mC).

---

## 3. Parte III — Implementação do classificador (Colab)

- [x] **3.1** Selecionar e justificar o dataset → GSE186575 (D4)
- [x] **EDA** em [notebooks/exploracao/01_eda_gse186575.ipynb](notebooks/exploracao/01_eda_gse186575.ipynb) — principais achados:
  - 98 amostras (CRC 32 · STAD 35 · NC 31) × 58.288 promotores (CPM); sem NaN; ~9% zeros
  - **Sexo** domina o PC2 (STAD 77% masculino) → remover chrX/Y/MT (D6) → 53.768 promotores autossômicos
  - **PC1 (46%) técnico**, correlacionado com fração de zeros (r = −0,63); 2 CRC `mix` são outliers extremos
  - **Lote `mix`**: RF separa NC mix vs NC não-mix com AUC ~0,82; 25/31 controles são mix vs 2/67 cânceres → detecção parcialmente confundida com lote (D8)
  - **Idade**: controles ~10 anos mais jovens → outro confundidor da detecção
- [ ] **3.2** Montar matriz amostras × features (promotores), definir rótulos; pré-processamento: log2(CPM+1), remoção de chrX/Y/MT, filtro de promotores de baixa contagem/variância
- [ ] **3.3** Seleção/redução de features **dentro do pipeline** (ajustada só no treino de cada fold) — ex.: filtro por variância (não supervisionado) + `SelectKBest`/teste t (supervisionado) dentro de `sklearn.Pipeline`
- [ ] **3.4** Treinar ≥1 algoritmo — **Random Forest como modelo principal** (decisão D3). Opcional: regressão logística simples só como *baseline* de comparação, citando que é a família de modelo usada no artigo.
  - RF: `class_weight="balanced"` (desbalanceamento), `n_estimators` ~500, ajustar `max_features`/`max_depth` via CV interna
  - Importância de features: *impurity importance* + *permutation importance* (mais confiável) → responde a pergunta 12
- [ ] **3.5** Validação: hold-out estratificado + CV estratificada no treino; `GroupKFold`/`StratifiedGroupKFold` **por indivíduo** se houver múltiplas amostras por pessoa (ex.: TCGA tumor + normal do mesmo paciente)
- [ ] **3.6** Métricas: matriz de confusão, precisão, sensibilidade (recall), especificidade, F1; acurácia e ROC-AUC quando apropriadas (macro/por classe no multiclasse)
- [ ] **(Extensão opcional)** Fluxo hierárquico: detecção (câncer vs normal) → TOO nos positivos, imitando o artigo

Estrutura sugerida do notebook:
1. Setup e download dos dados (reprodutível, URLs fixas)
2. Registro das fontes, classes e estratégia de validação
3. EDA (distribuição de classes, NAs, PCA/UMAP exploratório)
4. Split treino/teste (antes de qualquer seleção supervisionada)
5. Pipeline (pré-proc + seleção + modelo) e CV
6. Avaliação no teste
7. Features mais importantes (coeficientes / importância) e anotação de genes
8. Conclusões

---

### 3.7 Vazamento × confundimento (conceito para o relatório)
- **Vazamento (data leakage)** = informação do conjunto de teste influenciando o treino. Exemplos: selecionar features ou normalizar com o dataset inteiro antes do split; mesmo paciente no treino e no teste. **É erro de procedimento → é evitado 100%** com `Pipeline` + seleção dentro de cada fold + split por indivíduo.
- **Confundimento / efeito de lote (batch effect)** = variável técnica ou demográfica (lote `mix`, sexo, idade) correlacionada com o rótulo. O modelo pode aprender esse "atalho" (*shortcut learning*). **É característica dos dados, não erro de procedimento** → não dá para eliminar totalmente; dá para **reduzir, medir e discutir**.
- Por que não "padronizar" o lote (ex.: ComBat): quando o lote está quase perfeitamente alinhado com a classe (25/31 NC são mix; 0/35 STAD são mix), corrigir o lote **apaga junto o sinal biológico** — não dá para separar o que é lote do que é câncer. E, se feita no dataset inteiro antes do split, a correção vira **vazamento**. O que fazemos: remover o confundidor quando ele está nas features (sexo → tirar chrX/Y), e para o lote, análise de sensibilidade nas amostras que cruzam o lote (6 NC não-mix, 2 CRC mix).
- Nota: "mix" é uma **inferência** nossa a partir do nome das amostras; o GEO não documenta o que significa.

## 4. Parte IV — Discussão (responder no relatório)

- [ ] 11. Há diferenças de metilação que separam as classes?
- [ ] 12. Quais features/regiões mais contribuem? (genes associados)
- [ ] 13. Desempenho geral e por classe
- [ ] 14. Limitações: tamanho amostral, desbalanceamento, alta dimensionalidade (p ≫ n)
- [ ] 15. Diferenças vs metodologia original (material, tecnologia, escala, arquitetura, threshold de especificidade)
- [ ] 16. Por que não é um teste clínico de detecção precoce (dados de tecido/sintomáticos, sem coorte assintomática, sem validação prospectiva, sem PPV populacional, etc.)

---

## 5. Registro de decisões (justificativas)

> Cada decisão relevante entra aqui com data, o que foi decidido, alternativas consideradas e o porquê. O professor pede isso explicitamente — este log vira base para o relatório.

| # | Data | Decisão | Alternativas | Justificativa |
|---|---|---|---|---|
| D1 | 2026-11-17 | Manter este `.md` como plano vivo + log de decisões | — | Acompanhamento das etapas e rastreabilidade das escolhas para o relatório |
| D2 | 2026-11-17 | Não usar dados CCGA/STRIVE | Solicitar acesso à GRAIL | Sem accession público; STRIVE declara `ipdSharing: NO`; artigo sem *data availability*. Inviável no prazo da disciplina |
| D3 | 2026-11-17 | Random Forest como algoritmo principal | Regressão logística regularizada (usada no artigo), SVM | Domínio do autor para explicar/apresentar; lida bem com p ≫ n e não-linearidade; dá importância de features nativamente. Diferença em relação ao artigo será discutida na pergunta 15 |
| D4 | 2026-11-17 | **(aprovado)** GSE186575 como dataset principal | GSE81314, TCGA, GSE275598 | Único candidato público com **plasma cfDNA + controles saudáveis + >1 tipo de câncer** com n razoável → permite reproduzir o fluxo detecção → TOO do artigo |
| D5 | 2026-11-17 | *(opcional, decidir depois)* TCGA por coorte como experimento complementar | PanCan único (14 GB) | Arquivo PanCan grande demais para o Colab; por coorte cabe na memória |
| D6 | 2026-11-17 | Remover promotores de chrX, chrY, MT e sem anotação | Manter todos; ajustar por sexo | PCA mostrou que o sexo domina o PC2; STAD é 77% masculino → o modelo aprenderia sexo em vez de tecido de origem |
| D7 | 2026-11-17 | Versionar no git os dados do GSE186575 (22 MB) e a anotação Ensembl (3 MB); TCGA fica fora do git | Baixar tudo sempre do GEO | Garante que o notebook roda no dia da apresentação mesmo se o GEO estiver fora; o notebook ainda baixa do GEO se faltar algum arquivo. Arquivos do TCGA têm centenas de MB (limite do GitHub: 100 MB/arquivo) |
| D8 | 2026-11-17 | Tratar o lote `mix` como confundidor: reportar e fazer análise de sensibilidade (avaliar à parte os 6 NC não-mix e os 2 CRC mix) | Descartar amostras mix (sobrariam só 6 controles) | Remover os mix inviabiliza a tarefa de detecção; é melhor manter, medir e discutir o viés |
| D10 | 2026-11-17 | TCGA (tecido, 450K) como dataset principal; GSE186575 como tentativa documentada | Manter GSE186575 principal + TCGA adicional | Detecção em plasma ≈ acaso (AUC 0,55), sem promotores significativos após FDR e sinal de lote > sinal de câncer; TCGA permite experimento viável, validação agrupada por paciente (item 9 do enunciado) e fluxo hierárquico. O problema passa a ser explicitamente **classificação de tecido**, não detecção em cfDNA |
| D11 | 2026-11-17 | TCGA: até 150 tumores/coorte + todos os normais + tumores pareados; filtro não supervisionado (sem NaN, autossomos, top 20 mil por variância) fora da CV | Usar todas as amostras / seleção dentro da CV | Equilibra classes (BRCA tinha 790 tumores) e cabe no Colab gratuito; enunciado só exige ajuste no treino para seleção **supervisionada** |
| D9 | 2026-11-17 | Estrutura de pastas + notebooks gerados por script (`scripts/build_*.py`) | Editar o .ipynb à mão | Diff legível no git; notebook de entrega único em `notebooks/trabalho_final.ipynb`, rascunhos em `notebooks/exploracao/` |

---

## 6. Referências (ir acumulando)

- Liu, M. C. et al. (2020). Annals of Oncology, 31(6), 745–759. doi:10.1016/j.annonc.2020.02.011
- ClinicalTrials.gov: NCT02889978 (CCGA), NCT03085888 (STRIVE)
- GSE186575 — PMID 37992683
- GSE81314 — Li W. et al. (2017) *Cell Research* 27:1243–1257
- Moss J. et al. (2018) *Nat Commun* 9:5068 (GSE122126)
- Zeng et al. (2026) *Nature Cancer* — doi:10.1038/s43018-026-01116-3 (exemplo de acesso controlado)
- The Cancer Genome Atlas, via UCSC Xena (Goldman et al. 2020, *Nat Biotechnol*)

---

## Status atual

- [x] Leitura do enunciado e do artigo
- [x] Criação deste plano
- [x] Parte II — investigação dos datasets (tabela 2.2 preenchida)
- [x] Estrutura do repositório + git local + EDA do GSE186575
- [ ] Criar repositório no GitHub e ajustar `REPO_URL` nos notebooks
- [x] Notebook de modelagem (rascunho) [02_modelagem_rf.ipynb](notebooks/exploracao/02_modelagem_rf.ipynb) — resultados:
  - **Detecção:** nested CV AUC 0,55 ± 0,18 (≈ acaso); especificidade 0,32. Idade sozinha: AUC 0,70
  - **TOO CRC vs STAD:** nested CV AUC 0,79 ± 0,14; hold-out AUC 0,87, acurácia 0,81. Idade sozinha: AUC 0,54 (não explica)
  - **Hierárquico (3 classes):** acurácia 0,40
  - **Lote:** os 6 NC não-mix foram todos chamados de câncer
  - **Sinal univariado (FDR 5%):** 0 promotores em câncer vs saudável, 0 em CRC vs STAD; 526 em NC mix vs NC não-mix (só lote)
  - **Variações de pré-processamento** (protein-coding, rank, z-score por amostra): detecção no máximo AUC 0,64
  - **Demo de vazamento (rótulos aleatórios):** AUC 0,84 com seleção fora da CV vs 0,27 dentro → seleção fora da CV inventa sinal
  - Top features: muitos pseudogenes (RPL/RPS…P) → provável ruído técnico de mapeamento
- [x] **Decisão D10:** TCGA vira o dataset **principal**; GSE186575 fica documentado como tentativa em plasma (motivo: detecção ≈ acaso, sinal dominado por lote/idade)
- [x] Notebook TCGA (rascunho) [03_tcga_rf.ipynb](notebooks/exploracao/03_tcga_rf.ipynb) — 6 tecidos, 1.290 amostras (910 tumor / 380 normal), 919 pacientes (371 com par tumor+normal), 20 mil CpGs; nested CV agrupada por paciente:
  - **Etapa 1 tumor vs normal:** acurácia 0,984 · AUC 0,997 · sens 0,987 · espec 0,979 (estômago só tem 2 normais → espec. 0/2 nesse tecido)
  - **Etapa 2 TOO (6 classes):** acurácia 0,980 · F1 macro 0,980 · AUC macro 0,999; pior par: cólon ↔ estômago (9 cólon → estômago)
  - **Hierárquico (7 classes):** acurácia 0,971
- [ ] **PRÓXIMO:** montar `notebooks/trabalho_final.ipynb` (TCGA principal + seção curta do GSE186575) e preencher resultados no relatório
