# Plano de Execução — Detecção e classificação de câncer por metilação do DNA

**Disciplina:** Genômica Computacional
**Artigo de referência:** Liu, M. C. et al. (2020). *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*. Annals of Oncology, 31(6), 745–759. DOI: 10.1016/j.annonc.2020.02.011
**Documentos-fonte:** `docs/atividade_genomica_methylacao_multicancer.pdf` (enunciado) e `docs/artigo_referencia.pdf` (artigo)
**Última atualização:** 2026-11-18

> Legenda: `[ ]` pendente · `[~]` em andamento · `[x]` concluído
> Este arquivo é **autossuficiente**: dá para estudar o trabalho todo só por ele, sem abrir os notebooks.

---

## 0. Resumo executivo (onde estamos)

**O trabalho em uma frase:** classificar câncer e tecido de origem pela metilação do DNA, com um **experimento principal em tecido tumoral (TCGA)** e uma **tentativa documentada em plasma (GSE186575)** que serve de contraste com o artigo.

**A história em 6 frases:**
1. O artigo detecta câncer e o tecido de origem (TOO) a partir da metilação do cfDNA do plasma, mas usa dados **privados** da GRAIL (CCGA/STRIVE).
2. **Experimento principal — TCGA:** metilação de **tecido**, 6 tipos de câncer, 1.290 amostras de 919 pacientes. O Random Forest, validado por nested CV agrupada por paciente, separa tumor de normal e identifica o tecido de origem com acurácia de 0,98 (seção 4).
3. Uma **auditoria de vazamento** confirmou que o resultado é legítimo: a tarefa em tecido é genuinamente fácil, porque o sinal está espalhado pelo genoma.
4. **Tentativa documentada — plasma (GSE186575):** foi o primeiro dataset escolhido, por ser o mais parecido com o artigo (colorretal, gástrico e saudáveis). A **detecção falhou** (AUC 0,55 ≈ acaso), dominada por **lote técnico e idade** (seção 5).
5. Por isso o TCGA virou o experimento principal (decisão D10), e o plasma ficou como evidência e contraste.
6. O **contraste tecido × plasma** é o eixo da discussão. O sinal existe na fonte (o tumor), mas fica diluído no sangue. Por isso o artigo precisou de escala enorme e de sequenciamento direcionado.

**Números-chave:**

| Experimento | Tarefa | Resultado |
|---|---|---|
| **Principal — tecido (TCGA)** | Tumor vs normal | Acurácia **0,984** · AUC **0,997** |
| **Principal — tecido (TCGA)** | Tecido de origem (6 classes) | Acurácia **0,980** · F1 macro **0,980** |
| **Principal — tecido (TCGA)** | Hierárquico (normal + 6 tecidos) | Acurácia **0,971** |
| Tentativa — plasma (GSE186575) | Câncer vs saudável | AUC **0,55** ± 0,18 (acaso) |
| Tentativa — plasma (GSE186575) | Colorretal vs gástrico | AUC **0,79** ± 0,14 |

**Próximo passo:** validar o notebook final no Colab e reestruturar o relatório (ver seção 12).

---

## 1. Entregáveis

- [~] **Relatório técnico** (`relatorio/relatorio.md`, depois LaTeX/Overleaf → PDF). O rascunho existe, mas ainda trata o GSE186575 como principal: **precisa ser reestruturado** (ver 12).
- [x] **Notebook Colab reproduzível único** (`notebooks/trabalho_final.ipynb`, gerado por `scripts/build_final_notebook.py`): TCGA como experimento principal (com auditoria), seção de tentativa em plasma, comparação e conclusões. Executado de ponta a ponta localmente; falta testar no Colab.
- [x] **Registro de fontes, classes e validação:** seção 10 deste arquivo (vai para o relatório e para o início do notebook).

Critério central do professor: **escolha de dados adequada** e **prevenção de vazamento de informação**.

---

## 2. Parte I — Estudo do artigo (conteúdo validado no PDF)

### 2.1 Fundamentação biológica
- **Metilação do DNA:** adição de um grupo metil à citosina (5mC), quase sempre em dinucleotídeos **CpG**. O genoma tem cerca de 28 a 30 milhões de CpGs. **Ilhas CpG** em promotores costumam ficar desmetiladas, e a metilação de um promotor tende a **silenciar** o gene. O padrão é estável e **específico de cada tipo celular**, funcionando como uma "impressão digital" do tecido.
- **Câncer:** apresenta **hipometilação global** (instabilidade genômica) e **hipermetilação focal** de promotores (silencia genes supressores de tumor e de reparo). São alterações precoces e recorrentes, por isso bons biomarcadores.
- **cfDNA:** fragmentos de cerca de 167 pb no plasma, liberados por apoptose e necrose. Em saudáveis, vem principalmente de células do sangue. **ctDNA** é a fração de origem tumoral, muitas vezes **menor que 1%** em estádios iniciais. **Biópsia líquida** é analisar esse material por uma coleta de sangue.
- **Por que a metilação serve para detecção e para TOO:**
  - o ctDNA carrega padrões aberrantes que não aparecem em saudáveis (**detecção**);
  - esses padrões também refletem o tecido de origem (**TOO**).
  - Mutações não indicam o tecido e são confundidas por **CHIP** (hematopoiese clonal).

### 2.2 Amostras e dados do artigo

| Conjunto | Câncer | Não-câncer (CCGA + STRIVE) | Total |
|---|---|---|---|
| Treinamento | 1.742 | 2.978 (1.391 + 1.587) | 4.720 |
| Validação | 740 | 1.229 (614 + 615) | 1.969 |
| **Total** | **2.482** | **4.207** | **6.689** |

- 4.316 participantes foram analisáveis na análise primária (3.052 no treino e 1.264 na validação). 354 amostras do CCGA ficaram reservadas como referência de biópsias tumorais.
- **CCGA (NCT02889978):** prospectivo, multicêntrico, caso-controle, 15.254 participantes de 142 centros. Tem 3 sub-estudos: (1) descoberta, que comparou WGBS, WGS e mutações e elegeu a metilação; (2) **este artigo**, com treino e validação; (3) validação adicional. Forneceu os casos de câncer e parte dos controles.
- **STRIVE (NCT03085888):** 99.286 mulheres em mamografia de rastreio, em 35 centros. Forneceu **controles extras** (de um único centro) para treinar e medir especificidade acima de 99% com mais de 90% de confiança.
- **Tecnologia:** cfDNA (até 75 ng), conversão por **bissulfito**, captura híbrida com painel customizado (Twist) e sequenciamento NovaSeq 150 pb *paired-end*. Mediana de 113 M fragmentos (139× único). Lê o estado de metilação de **cada CpG em cada fragmento**.
- **Regiões:** escolhidas com WGBS do 1º sub-estudo e **arrays de metilação do TCGA**. São **103.456 regiões (17,2 Mb), cerca de 1,1 milhão de CpGs**.
- **Separação e vazamento:**
  - Salvaguardas: validação independente; classificador **travado** antes de liberar a validação; análise **duplo-cega**; equipe de integridade de dados separada; distribuição pré-especificada por centro e estádio; análise de centros removidos; taxa de falsos positivos (FPR) do CCGA (0,7%) parecida com a do STRIVE (0,6%).
  - Riscos residuais:
    - STRIVE = só mulheres e de um único centro (confundidor possível);
    - regiões escolhidas com dados de descoberta;
    - **várias amostras por participante** (11.154 amostras de 5.854 pessoas), o que exige separar por indivíduo;
    - acompanhamento incompleto dos controles.

### 2.3 Estratégia de classificação

| | (i) Detecção | (ii) Tecido de origem |
|---|---|---|
| Entrada | Todas as amostras | Só as positivas em (i) |
| Saída | Câncer / não-câncer | Um entre cerca de 20 sítios |
| Métrica | Especificidade (~99%) + sensibilidade | Acurácia do TOO |

- **Arquitetura:** *source models* por região (o quanto o padrão se parece com cada tipo de câncer) → features → **par de regressões logísticas em ensemble** (uma para detecção, outra para TOO). Limiar fixado para especificidade de cerca de 99%.
- **Resultados (validação):**
  - especificidade de **99,3%**;
  - sensibilidade nos estádios I–III de **67,3%** (12 cânceres pré-especificados) e **43,9%** (todos);
  - sensibilidade cresce com o estádio: 18% (I), 43% (II), 81% (III), 93% (IV);
  - TOO previsto em 96% dos positivos, com **93%** de acurácia.

---

## 3. Parte II — Datasets públicos

### 3.1 CCGA/STRIVE: não disponíveis
- O artigo **não tem** declaração de disponibilidade de dados nem número de acesso (GEO, EGA ou dbGaP).
- ClinicalTrials.gov: **STRIVE `ipdSharing: NO`**; o CCGA não declara plano de compartilhamento. Os dados são proprietários da GRAIL.
- Exemplo de **acesso controlado:** Zeng et al. (2026, *Nature Cancer*), com 1.294 amostras de cfMeDIP-seq. Os dados estão no EGA e exigem um *Data Access Agreement* com comitê institucional, inviável no prazo.

### 3.2 Candidatos (verificados nas APIs do GEO e do Xena em 17/11/2026)

| Campo | **GSE186575** | GSE81314 | **TCGA (Xena)** | GSE314062 | GSE275598 | GSE122126 |
|---|---|---|---|---|---|---|
| Publicação | PMID 37992683 | Li et al. 2017 *Cell Res* | TCGA/PanCanAtlas | PMID 42524278 | PMID 39187719 | Moss et al. 2018 |
| Amostras | 98 | 64 pacientes úteis | ~100–900/coorte | 5 | 93 | 59 cfDNA |
| Classes | CRC 32 · gástrico 35 · saudável 31 | 7 cânceres (4–15 cada) + 8 saudáveis + 7 HBV | Tumor vs normal por tipo | 3 pacientes HNC | Tumores cerebrais pediátricos | Sepse, transplante, poucos cânceres |
| Material | **Plasma** | Plasma | **Tecido** | Plasma | Plasma | Plasma + células |
| Tecnologia | cfMeDIP-seq | 5hmC (não 5mC) | 450K (bissulfito) | Bissulfito direcionado | cfMeDIP-seq | 450K/EPIC |
| Formato | CPM por promotor | Por gene | Beta por CpG | Bruto | CPM | IDATs |
| Viabilidade | Alta (tentado) | Média | **Alta (escolhido)** | Inviável | Fora do tema | Inviável |

- Os 4 sugeridos pelo professor foram TCGA, GSE186575, GSE314062 e GSE81314. O GSE275598 e o GSE122126 vieram da nossa busca.

### 3.3 TCGA: amostras disponíveis por coorte (450K)

| Coorte | Tumor | Normal | Tamanho |
|---|---|---|---|
| BRCA (mama) | 790 | 98 | 783 MB |
| COAD (cólon) | 299 | 38 | 309 MB |
| STAD (estômago) | 396 | **2** | 360 MB |
| LUAD (pulmão) | 460 | 32 | 450 MB |
| LIHC (fígado) | 379 | 50 | 386 MB |
| KIRC (rim) | 320 | 160 | 423 MB |

Outras coortes conferidas: LUSC, READ, PAAD, PRAD, HNSC, THCA, UCEC e ESCA. O arquivo PanCan único tem **14 GB**, inviável no Colab gratuito.

---

## 4. Experimento principal — Tecido tumoral (TCGA)

### 4.1 Dados e amostragem
- **Fonte:** UCSC Xena, arquivos `TCGA.<COORTE>.sampleMap/HumanMethylation450.gz`, mais o *probe map* `illuminaMethyl450_hg19_GPL16304_TCGAlegacy` (CpG → gene/cromossomo).
- **Medida:** **beta value** (0 a 1) por CpG, ou seja, a fração das cópias de DNA metiladas naquele CpG. Química por bissulfito, como no artigo.
- **6 tecidos:** cólon (COAD), estômago (STAD), pulmão (LUAD), fígado (LIHC), mama (BRCA) e rim (KIRC).
- **Código TCGA:** `TCGA-XX-YYYY-TT` → paciente = `TCGA-XX-YYYY`; `01` = tumor primário, `11` = normal adjacente.
- **Amostragem (D11):** todos os normais, os tumores **pareados** desses normais e tumores aleatórios até 150 por coorte (semente 42).

| Tecido | Normal | Tumor | Total |
|---|---|---|---|
| cólon | 38 | 150 | 188 |
| estômago | 2 | 150 | 152 |
| fígado | 50 | 150 | 200 |
| mama | 98 | 150 | 248 |
| pulmão | 32 | 150 | 182 |
| rim | 160 | 160 | 320 |
| **Total** | **380** | **910** | **1.290** |

São **919 pacientes**, dos quais **371** têm tumor e normal pareados.

### 4.2 Pré-processamento (`notebooks/trabalho_final.ipynb`, seção 3; versão de testes em `notebooks/exploracao/03_tcga_rf.ipynb`)
1. Mantidos só os CpGs autossômicos: 385.809 de 395.985 (remove o sinal de sexo).
2. Juntadas as 6 coortes; removidos os 49.549 CpGs com algum NaN, restando 336.260.
3. **Filtro não supervisionado** (não usa rótulos): os **20 mil CpGs de maior variância**. O enunciado só exige ajuste no treino para seleção **supervisionada**, e a auditoria mostrou que fazer esse filtro dentro da CV não muda nada.
4. Cache da matriz completa (1.290 × 336.260, sem NaN) em `data/processed/tcga_6coortes_completa.parquet` (fora do git). O filtro de variância é refeito a partir dela, e a auditoria usa a mesma matriz.

**PCA:** PC1 explica 34,7% e PC2 11,6%. Os **normais formam grupos compactos por tecido** (rim, mama, pulmão, fígado bem separados). Os **tumores se espalham**, refletindo heterogeneidade e desregulação da metilação. O fígado tumoral forma um "braço" próprio.

### 4.3 Método
- **Pipeline:** `SelectKBest(f_classif)` **dentro** do pipeline (supervisionado), seguido de RF com 500 árvores e `class_weight="balanced"`.
- **Validação:** **nested CV com `StratifiedGroupKFold`, agrupada por paciente**.
  - Externa com 5 folds; interna com 3 folds escolhendo k ∈ {100, 500, 2000}. O k = 2000 foi escolhido na maioria dos folds.
  - Tumor e normal do mesmo paciente **nunca** ficam um no treino e outro no teste. Isso cumpre o item 9 do enunciado e evita que o modelo reconheça o indivíduo em vez do tumor.
- **Etapa 1:** tumor vs normal. **Etapa 2:** tecido de origem entre os tumores (6 classes). **Hierárquico:** etapa 1 e, nos positivos, a etapa 2.
- **Importância das features:** *permutation importance* num split agrupado 75/25 (k = 500), anotada com o gene pelo *probe map*.

### 4.4 Resultados

**Etapa 1 — tumor vs normal:** acurácia **0,984** · ROC-AUC **0,997**

| Classe | n | Precisão | Sensibilidade | Especificidade | F1 |
|---|---|---|---|---|---|
| normal | 380 | 0,969 | 0,979 | 0,987 | 0,974 |
| tumor | 910 | 0,991 | 0,987 | 0,979 | 0,989 |

Por tecido:

| Tecido | Sensibilidade | Especificidade |
|---|---|---|
| cólon | 1,000 | 0,921 |
| estômago | 0,993 | **0,000 (só 2 normais)** |
| fígado | 0,967 | 0,960 |
| mama | 0,993 | 0,990 |
| pulmão | 0,987 | 1,000 |
| rim | 0,981 | 1,000 |

**Etapa 2 — tecido de origem:** acurácia **0,980** · F1 macro **0,980** · AUC macro (OvR) **0,999**

| Tecido | Precisão | Sensibilidade | Especificidade | F1 |
|---|---|---|---|---|
| cólon | 0,979 | 0,940 | 0,996 | 0,959 |
| estômago | 0,929 | 0,960 | 0,986 | 0,944 |
| fígado | 1,000 | 1,000 | 1,000 | 1,000 |
| mama | 0,987 | 1,000 | 0,997 | 0,993 |
| pulmão | 0,993 | 0,987 | 0,999 | 0,990 |
| rim | 0,994 | 0,994 | 0,999 | 0,994 |

**Hierárquico (7 classes, out-of-fold):** acurácia **0,971**
- 372 de 380 normais corretos. Os 8 normais chamados de tumor foram classificados como cólon (3), estômago (2), fígado (2) e mama (1).
- 12 tumores chamados de normal: fígado 5, rim 3, pulmão 2, estômago 1, mama 1.
- **Maior confusão: 9 tumores de cólon chamados de estômago** (e 3 de estômago chamados de cólon). Os dois são órgãos do trato gastrointestinal.

**CpGs mais relevantes (etapa 2), beta médio por tecido:**
- **KANK1** (2 CpGs): metilado no cólon (0,67–0,74), intermediário no estômago (0,23–0,27) e desmetilado nos demais (≤0,16). É um **marcador de cólon**.
- **ERI3**: desmetilado no cólon (0,09), intermediário no estômago (0,46) e alto nos demais (0,63–0,87).
- **RBP5, LITAF, RP1-60O19.1:** **desmetilados no fígado** (0,31–0,39) e altos nos demais. São **marcadores de fígado**. O RBP5 também é baixo no rim.
- **TCF3, MYBPC3, cg10900049:** baixos no cólon (0,25–0,44) e altos nos demais.
- Padrão geral: **o estômago fica "entre" o cólon e os outros tecidos**, o que explica a confusão cólon ↔ estômago.
- Etapa 1 (tumor vs normal): B3GNTL1, FAM78B, HMGA2, FMNL2, MIR646HG, COL8A1 e outros.
- Os valores de *permutation importance* são baixos (≤0,007) porque há **muitos CpGs redundantes**: tirar um só não derruba o modelo. Isso é coerente com a auditoria (o sinal está espalhado).
- Tabelas completas em `results/tables/tcga_importancias_etapa{1,2}.csv`.

### 4.5 Auditoria de vazamento (`notebooks/trabalho_final.ipynb`, seção 3.8; versão em script: `scripts/auditoria_tcga.py`)

| Teste | Etapa 1 | Etapa 2 | Conclusão |
|---|---|---|---|
| Referência | acc 0,982 · AUC 0,997 | acc 0,978 · AUC 0,999 | — |
| IDs duplicados / pacientes em treino e teste | 0 / 0 | 0 / 0 | Agrupamento correto |
| Filtro de variância **só no treino** | acc 0,982 | acc 0,977 | O filtro global não inflou nada |
| **20 mil CpGs aleatórios** | acc 0,988 | acc 0,980 | O sinal está espalhado pelo genoma |
| **Rótulos embaralhados** | AUC 0,504 | AUC 0,503 (acc 0,18 ≈ 1/6) | O pipeline não inventa sinal |
| Árvore com 3 níveis | acc 0,970 | acc 0,647 | A etapa 1 é trivial; a etapa 2 precisa de muitos CpGs |
| **1 único CpG** | acc 0,905 | acc 0,331 | Um CpG separa 90% de tumor vs normal |

**Risco que não dá para testar:** **lote ou placa** de processamento do TCGA confundido com o tipo de câncer, porque os IDs do Xena não trazem a placa. É mitigado pela coerência biológica (normais agrupados por tecido, marcadores tecido-específicos) e entra como **limitação**.

---

## 5. Tentativa documentada — Plasma (GSE186575)

### 5.1 Dados
- cfMeDIP-seq (anticorpo que "pesca" o DNA metilado) de plasma, Peking University First Hospital.
- Matriz: **98 amostras × 58.288 promotores**, em CPM (fragmentos metilados por milhão). Sem NaN; cerca de 9,3% de zeros.
- **Feature** = quantidade de DNA metilado no promotor de um gene. É sinal agregado, não % por CpG.
- Uma amostra por pessoa.

### 5.2 Análise exploratória (`notebooks/exploracao/01_eda_gse186575.ipynb`; resumo no notebook final, seção 4.1)
- **Sexo dominava o PC2:** 44 de 58 homens ficam no cluster "PC2 alto" e 33 de 39 mulheres no "baixo". O grupo gástrico é 77% masculino, então o modelo poderia separar gástrico de colorretal pelo sexo. **Removemos os promotores de chrX, chrY, MT e sem anotação**, restando 53.768. Depois disso, o efeito some.
- **PC1 (46%) é técnico:** se correlaciona com a fração de zeros (r = −0,63). Os 2 CRC "mix" são outliers extremos.
- **Lote "mix":**
  - 25 de 31 controles são "mix", contra 2 de 32 CRC e 0 de 35 gástricos;
  - um RF distingue controles mix de não-mix (todos saudáveis) com AUC de cerca de 0,82;
  - o significado de "mix" **não está documentado no GEO**; é uma hipótese nossa.
- **Idade:** os controles são cerca de 10 anos mais jovens (mediana de 50 anos, contra 62–64).

### 5.3 Método (`notebooks/exploracao/02_modelagem_rf.ipynb`; no notebook final, seção 4.2)
- **Pré-processamento:** `log2(CPM+1)`, amostra a amostra, mais a remoção de X/Y/MT.
- **Pipeline:**
  1. filtro de presença (sinal > 0 em ≥50% das amostras de treino);
  2. `SelectKBest(f_classif)`;
  3. RF com 500 árvores e `class_weight="balanced"`.
- **Validação:** nested CV (externa com 5 folds × 3 repetições; interna com 5 folds escolhendo k ∈ {100, 500, 2000} e `max_features` ∈ {sqrt, 0,1}), mais um hold-out 70/30.

### 5.4 Resultados

| Tarefa | Nested CV (média ± desvio) | Hold-out |
|---|---|---|
| Detecção (67 câncer vs 31 saudável) | AUC 0,553 ± 0,182 · acc 0,623 · sens 0,763 · **espec 0,324** | AUC 0,595 · acc 0,667 |
| TOO (32 CRC vs 35 STAD) | AUC 0,789 ± 0,136 · acc 0,711 · F1 0,726 | AUC 0,873 · acc 0,810 |
| Hierárquico (3 classes, out-of-fold) | acurácia 0,40 | — |

**Investigação da falha:**
- **Idade sozinha:** AUC de **0,70** na detecção (melhor que a metilação!) e 0,54 no TOO (não explica o TOO).
- **Lote:** os **6 controles não-mix foram todos chamados de câncer**.
- **Teste por promotor com FDR de 5%:**
  - 0 promotores significativos em câncer vs saudável;
  - 0 em CRC vs gástrico;
  - **526 em controle mix vs controle não-mix** → o lote deixa mais marca do que o câncer.
- **Outros pré-processamentos** (só genes codificadores, rank e z-score por amostra): a detecção chega no máximo a AUC de 0,64.
- **Features mais importantes:** muitos **pseudogenes** (RPL/RPS…P), sinal de ruído técnico de mapeamento.

### 5.5 Por que desistimos (decisão D10)
A detecção ficou no acaso, não havia nenhum promotor significativo e o sinal de lote era maior que o de câncer. Isso **não permite um experimento de classificação viável**, que é o mínimo do enunciado. A tentativa fica no trabalho como **evidência** de que o problema em plasma é difícil e que exige a escala e a tecnologia do artigo.

### 5.6 Demonstração de vazamento (notebook final, seção 4.3)
Com **rótulos aleatórios**:
- seleção de features feita **antes** da CV: AUC **0,84** (sinal inventado);
- seleção **dentro** do pipeline: AUC **0,27**.

O valor abaixo de 0,5 é flutuação esperada com n pequeno e rótulos sem sinal; o ponto é que não há sinal espúrio.

---

## 6. Comparação plasma × tecido (eixo da discussão)

| | Tecido (TCGA) | Plasma (GSE186575) | Artigo |
|---|---|---|---|
| Material | Tumor | cfDNA | cfDNA |
| Fração tumoral | Alta (maior parte do DNA) | < 1% em geral | < 1% em geral |
| Tecnologia | Bissulfito, array | Imunoprecipitação | Bissulfito, sequenciamento |
| Resolução | CpG (~485 mil) | Promotor (~58 mil) | CpG por fragmento (~1,1 milhão) |
| n | 1.290 | 98 | 6.689 |
| Detecção | AUC 0,997 | AUC 0,55 | Espec. 99,3% / sens. 55% (I–IV) |
| TOO | acc 0,98 (6 classes) | AUC 0,79 (2 classes) | acc 93% (~20 sítios) |

**Mensagem:** a assinatura de metilação dos tecidos é forte na fonte. No plasma, ela fica diluída no DNA do sangue e é facilmente superada por ruído técnico (lote, idade). O artigo vence essa distância com escala, profundidade e leitura molécula a molécula.

---

## 7. Parte IV — Respostas (rascunho para o relatório)

**11. Há diferenças de metilação que separam as classes?**
- **Em tecido, sim, e muito.** 20 mil CpGs aleatórios bastam, e um único CpG separa 90% de tumor vs normal.
- **Em plasma (dados públicos), não de forma detectável para câncer vs saudável:** nenhum promotor passou na correção para testes múltiplos. Para colorretal vs gástrico, há sinal multivariado fraco (AUC de 0,79).

**12. Quais regiões contribuem mais?**
- **No tecido:** CpGs **tecido-específicos**, como KANK1 (metilado no cólon), RBP5/LITAF (desmetilados no fígado) e ERI3/TCF3 (desmetilados no cólon). O estômago tem valores intermediários entre cólon e os demais.
- **No plasma:** majoritariamente pseudogenes, o que indica ruído técnico.
- A importância de cada CpG é baixa por causa da redundância.

**13. Desempenho e variação entre classes.**
- **No tecido:** fígado, mama, pulmão e rim ficam ≥0,99. O cólon é o pior (sensibilidade de 0,94), confundido com o estômago, que é o tecido vizinho. A especificidade do estômago na etapa 1 não pode ser avaliada (só 2 normais).
- **No plasma:** a detecção falha (especificidade de 0,32), e o TOO é moderado.

**14. Limitações (tamanho amostral, desbalanceamento, alta dimensionalidade).**
- **p ≫ n:** 20 mil a 53 mil features contra 98 a 1.290 amostras. É mitigado por seleção dentro da CV e nested CV.
- **Desbalanceamento:** estômago com 2 normais; tumores limitados a 150 por coorte; `class_weight` balanceado; métricas por classe.
- **Plasma:** n = 98, um único centro, lote e idade confundidos com a classe.
- **TCGA:** possível lote por coorte; os tumores são casos já diagnosticados (sem estádio controlado).
- Com 31 controles, é impossível medir uma especificidade de 99%.

**15. Diferenças em relação ao artigo.** Material (tecido ou MeDIP, contra bissulfito em cfDNA), escala, número de tipos, resolução, features (seleção por teste F contra *source models*), classificador (RF contra regressão logística), limiar padrão de 0,5 contra especificidade fixa de 99%, e validação cruzada contra conjunto independente, cego e travado.

**16. Por que não é um teste clínico.**
- Tecido não é sangue; o teste exigiria uma biópsia.
- O desenho é caso-controle, com pacientes já diagnosticados, e não uma população assintomática.
- Não há validação externa, prospectiva ou cega.
- A especificidade não foi medida na escala necessária, e não há valor preditivo positivo (PPV) populacional.
- Poucos tipos de câncer; sem doenças benignas como controle.
- Sem estratificação por estádio (a detecção precoce depende do estádio I).
- Sem demonstração de benefício clínico.

---

## 8. Defesa: perguntas prováveis e respostas

**"O resultado do TCGA não está bom demais? Teve vazamento?"**
Não. A auditoria (4.5) mostra: 0 pacientes cruzando folds, filtro só no treino igual, CpGs aleatórios igual, rótulos embaralhados dando AUC de 0,50. A tarefa é fácil porque o tecido tumoral é quase todo DNA tumoral e cada órgão tem um programa de metilação próprio (Hoadley et al., 2018).

**"Se é fácil, por que fazer?"**
1. É fácil para o classificador, não para um humano. Ninguém lê 485 mil valores.
2. Tem uso clínico real: **câncer de origem primária desconhecida** (CUP, 3–5% dos cânceres). Classificadores de metilação em tecido como o EPICUP (Moran et al., 2016) e o de tumores cerebrais de Capper et al. (2018) são usados para isso.
3. É o **controle positivo**: prova que o pipeline funciona quando o sinal existe. Com isso, a falha no plasma passa a ser atribuída aos dados, não ao método.
4. O próprio artigo usou arrays do TCGA para escolher as regiões do painel.
5. Nem tudo é fácil: cólon ↔ estômago, e a árvore com 3 perguntas só acerta 65% do TOO.

**"Por que trocar de dataset em vez de insistir no plasma?"**
O resultado não tinha sinal (0 promotores com FDR, AUC no acaso), e o lote dominava. Um experimento viável é o mínimo exigido. A tentativa foi mantida e documentada como evidência.

**"Por que Random Forest e não regressão logística como no artigo?"**
Lida bem com p ≫ n e interações, é pouco sensível à escala, dá importância de features e é interpretável na apresentação. A diferença de modelo é discutida na pergunta 15.

**"Por que o filtro de variância fora da CV?"**
É não supervisionado (não usa rótulos), e o enunciado só exige ajuste no treino para seleção supervisionada. Mesmo assim, testamos dentro da CV e o resultado foi idêntico.

**"Por que agrupar por paciente?"**
371 pacientes têm tumor e normal. Se um ficasse no treino e outro no teste, o modelo poderia reconhecer a metilação constitucional do indivíduo, e não o tumor.

**"O que é 'mix' no GSE186575?"**
Não está documentado. Inferimos um lote pelo nome e confirmamos estatisticamente (526 promotores diferem entre controles mix e não-mix).

**"Por que tirar os cromossomos X e Y?"**
No plasma, o sexo dominava o PCA, e o grupo gástrico era 77% masculino: o modelo aprenderia sexo em vez de tecido. No TCGA, mantivemos só os autossomos pelo mesmo motivo (mama é praticamente só feminina).

---

## 9. Registro de decisões

| # | Decisão | Alternativas | Justificativa | Status |
|---|---|---|---|---|
| D1 | Este `.md` como plano vivo + log de decisões | — | Rastreabilidade para o relatório | Vigente |
| D2 | Não usar dados CCGA/STRIVE | Pedir acesso à GRAIL | Sem número de acesso público; STRIVE `ipdSharing: NO`; artigo sem declaração de dados | Vigente |
| D3 | Random Forest como algoritmo | Regressão logística (artigo), SVM | Domínio do autor; p ≫ n; importância de features | Vigente |
| D4 | GSE186575 como dataset principal | GSE81314, TCGA | Único público com plasma + controles + 2 cânceres | **Substituída pela D10** |
| D5 | TCGA por coorte (não PanCan) | PanCan 14 GB | Cabe no Colab gratuito | Vigente |
| D6 | Remover chrX/Y/MT/sem anotação no plasma | Manter tudo | O sexo dominava o PC2; gástrico 77% masculino | Vigente |
| D7 | Versionar GSE186575 (22 MB) e Ensembl (3 MB) no git; TCGA fora | Baixar sempre | Robustez no dia da apresentação; limite de 100 MB/arquivo do GitHub | Vigente |
| D8 | Manter o lote "mix" e medir (análise de sensibilidade) | Descartar os mix (sobrariam 6 controles) | Descartar inviabilizaria a detecção | Vigente (só na tentativa com plasma) |
| D9 | Notebooks gerados por script (`scripts/build_*.py`) | Editar o .ipynb à mão | Diff legível no git | Vigente |
| D10 | **TCGA como principal**; GSE186575 como tentativa documentada | Plasma principal + TCGA extra | Detecção no plasma ≈ acaso; 0 promotores com FDR; lote > câncer. O TCGA permite um experimento viável, agrupamento por paciente e fluxo hierárquico. O problema passa a ser **classificação de tecido** | Vigente |
| D11 | TCGA: ≤150 tumores/coorte + todos os normais + pareados; filtro sem NaN + autossomos + top 20 mil por variância fora da CV | Todas as amostras; filtro dentro da CV | Equilibra classes (BRCA tinha 790 tumores); cabe no Colab; não supervisionado | Vigente (auditado) |
| D12 | Não aplicar correção de lote (ComBat) no plasma | ComBat | Lote quase alinhado à classe: a correção removeria o sinal biológico | Vigente |
| D13 | Manter a tentativa em plasma e a demonstração de vazamento no trabalho final | Omitir | Evidência para as perguntas 11, 14 e 16; mostra a prevenção de vazamento, critério central | Vigente |
| D14 | Relatório em LaTeX (Overleaf) → PDF | Word | Preferência do autor | Vigente |
| D15 | Não contrastar "vazamento × confundimento" como seção no relatório | Seção didática | O professor já domina a diferença; basta citar os confundidores | Vigente |

---

## 10. Registro de fontes, classes e validação (para o relatório e o notebook)

| Item | Experimento principal — TCGA (tecido) | Tentativa — GSE186575 (plasma) |
|---|---|---|
| Fonte | UCSC Xena: `TCGA.{COAD,STAD,LUAD,LIHC,BRCA,KIRC}.sampleMap/HumanMethylation450.gz` + probe map `illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | GEO FTP: `GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz` + `GSE186575_series_matrix.txt.gz`; anotação Ensembl BioMart (GRCh38) |
| Classes | Etapa 1: tumor (`01`) vs normal (`11`). Etapa 2: cólon, estômago, pulmão, fígado, mama, rim | Detecção: CRC+STAD vs NC. TOO: CRC vs STAD |
| Validação | Nested CV, `StratifiedGroupKFold` por paciente (externa 5 / interna 3) | Nested CV estratificada (externa 5×3 / interna 5) + hold-out 70/30; uma amostra por pessoa |
| Anti-vazamento | Seleção supervisionada no pipeline; agrupamento por paciente; auditoria | Seleção supervisionada no pipeline; demonstração com rótulos aleatórios |

---

## 11. Estrutura do repositório e reprodução

```
GenomicaComputacional/
├── PLANO_EXECUCAO.md            # este arquivo
├── README.md                    # visão geral + como rodar
├── requirements.txt             # numpy, pandas, scikit-learn, matplotlib, seaborn, pyarrow
├── docs/                        # enunciado e artigo (PDF)
├── data/raw/                    # GSE186575 (no git); data/raw/tcga/ (fora do git, ~2,7 GB)
├── data/reference/              # ensembl_genes_hg38.tsv
├── data/processed/              # cache parquet do TCGA, matriz completa (fora do git)
├── notebooks/exploracao/        # 01_eda_gse186575 · 02_modelagem_rf · 03_tcga_rf
├── notebooks/trabalho_final.ipynb  # ENTREGA: notebook único (gerado por scripts/build_final_notebook.py)
├── results/figures, tables/     # figuras e tabelas geradas
├── relatorio/relatorio.md       # rascunho do relatório
└── scripts/                     # build_final_notebook.py (entrega), build_*_notebook.py (testes), auditoria_tcga.py
```

- **Ambiente local:** `.venv` com Python 3.13. Tempo de execução: EDA do plasma ~1 min; modelagem do plasma ~12 min; TCGA ~7 min com os downloads prontos; auditoria ~10 min.
- **Git:** só local. Repositório no GitHub adiado; o envio pode ser feito de outras formas.
- **Download robusto:** o notebook do TCGA baixa para um arquivo `.part` e só renomeia no fim. Isso foi motivado pela queda de conexão durante o download do KIRC.

---

## 12. Status e próximos passos

**Concluído:**
- [x] Leitura do enunciado e do artigo; plano e log de decisões
- [x] Parte II: investigação de 6 datasets + CCGA/STRIVE
- [x] Experimento principal (TCGA): download, amostragem, nested CV agrupada, hierárquico, importâncias, heatmap
- [x] Tentativa com plasma (GSE186575): EDA, modelagem, investigação da falha, decisão de troca
- [x] Auditoria de vazamento do TCGA
- [x] Rascunho do relatório (Parte I, Parte II e metodologia do plasma)

**Pendente (em ordem):**
1. [x] **Notebook final** `notebooks/trabalho_final.ipynb`: registro de fontes/classes/validação → contexto e Parte II → **TCGA** (dados, PCA, pipeline, etapa 1, etapa 2, hierárquico, importâncias, auditoria) → **tentativa com plasma** (confundidores, FDR, nested CV, lote/idade, demonstração de vazamento) → comparação → conclusões.
2. [ ] **Validar o notebook final no Google Colab gratuito** (memória e tempo; o download do TCGA leva alguns minutos).
3. [ ] **Reestruturar o relatório:** TCGA como Parte III principal; plasma como subseção da Parte II ("avaliação experimental"); preencher resultados, Parte IV (seção 7), conclusão e as novas referências (Hoadley 2018, Moran 2016, Capper 2018).
4. [ ] Converter o relatório para LaTeX (Overleaf).
5. [ ] Montar os slides (seção 14).
6. [ ] (Opcional) Repositório no GitHub e ajuste de `REPO_URL` nos notebooks.

**Conferência deste plano (18/11/2026):**
- Todos os números das seções 4, 5 e 6 foram conferidos contra as saídas dos notebooks `02_modelagem_rf`, `03_tcga_rf`, as tabelas em `results/tables/` e a saída da `auditoria_tcga.py`.
- Os números do artigo (seção 2) foram conferidos no PDF.
- **Notebook final executado de ponta a ponta localmente (18/11/2026), sem erros.** Ele reproduz exatamente os números das seções 4 e 5:
  - TCGA: 0,984 / 0,997 · 0,980 / 0,980 / 0,999 · hierárquico 0,971;
  - auditoria idêntica;
  - plasma: AUC 0,553 / 0,789; FDR 0 / 0 / 526; idade 0,703; vazamento 0,838 × 0,265.
  - O notebook também mostra a tabela por classe do fluxo hierárquico: estômago tem F1 0,938, o pior; os demais ficam ≥ 0,95.

---

## 13. Referências

1. Liu, M. C. et al. (2020). *Annals of Oncology*, 31(6), 745–759. doi:10.1016/j.annonc.2020.02.011
2. ClinicalTrials.gov: NCT02889978 (CCGA); NCT03085888 (STRIVE)
3. GEO GSE186575 — PMID 37992683
4. Li, W. et al. (2017). *Cell Research*, 27, 1243–1257 (GSE81314)
5. Moss, J. et al. (2018). *Nature Communications*, 9, 5068 (GSE122126)
6. Zeng et al. (2026). *Nature Cancer*. doi:10.1038/s43018-026-01116-3
7. Goldman, M. J. et al. (2020). UCSC Xena. *Nature Biotechnology*, 38, 675–678
8. Hoadley, K. A. et al. (2018). Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell*, 173(2), 291–304
9. Moran, S. et al. (2016). Epigenetic profiling to classify cancer of unknown primary: a multicentre, retrospective analysis. *Lancet Oncology*, 17(10), 1386–1395
10. Capper, D. et al. (2018). DNA methylation-based classification of central nervous system tumours. *Nature*, 555, 469–474
11. Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32
12. Varma, S. & Simon, R. (2006). *BMC Bioinformatics*, 7, 91 (nested CV)
13. Ambroise, C. & McLachlan, G. J. (2002). *PNAS*, 99(10), 6562–6566 (viés de seleção)
14. Leek, J. T. et al. (2010). *Nature Reviews Genetics*, 11, 733–739 (efeito de lote)
15. Pedregosa, F. et al. (2011). Scikit-learn. *JMLR*, 12, 2825–2830

---

## 14. Protótipo da sequência de slides (apresentação do relatório técnico)

> Cerca de 15 a 18 slides, para uns 15 minutos. Entre colchetes, a figura ou tabela já gerada em `results/figures/` ou `results/tables/`.

| # | Título | Conteúdo | Visual |
|---|---|---|---|
| 1 | **Detecção de câncer por metilação do DNA** | Título, autor, disciplina, artigo de referência | — |
| 2 | **O problema** | Câncer é diagnosticado tarde; rastreamento existe para poucos tipos; proposta: um exame de sangue para vários cânceres e seu órgão de origem | Ícones: sangue → câncer? → órgão? |
| 3 | **Metilação do DNA e câncer** | CpG, ilhas CpG, silenciamento; hipometilação global e hipermetilação de promotores; padrão específico de cada tecido | Esquema de um CpG metilado/não metilado |
| 4 | **cfDNA e biópsia líquida** | Fragmentos no plasma; ctDNA < 1% em estádio inicial; por que a metilação dá detecção e tecido de origem | Esquema: tecidos → sangue → fragmentos |
| 5 | **O artigo: desenho do estudo** | CCGA + STRIVE; 6.689 participantes; tabela treino/validação; bissulfito + painel de 1,1 milhão de CpGs | Tabela 2.2 |
| 6 | **O artigo: classificação e resultados** | Duas etapas (detecção → TOO); *source models* + regressões logísticas; especificidade 99,3%; sensibilidade por estádio; TOO 93% | Tabela 2.3 + barras de sensibilidade por estádio |
| 7 | **Salvaguardas contra vazamento no artigo** | Validação travada e cega; equipe independente; centros removidos; CCGA × STRIVE; riscos residuais | Lista |
| 8 | **Busca de dados públicos** | CCGA/STRIVE indisponíveis (`ipdSharing: NO`); público × controlado; 6 candidatos avaliados; **escolha: TCGA (tecido) como principal, GSE186575 (plasma) como tentativa** | Tabela 3.2 resumida |
| 9 | **Experimento principal: TCGA** | 6 tecidos, 1.290 amostras, 919 pacientes, beta values 450K; tecido ≠ plasma (problema redefinido para classificação de tecido) | Tabela 4.1 |
| 10 | **Pipeline sem vazamento** | Filtro não supervisionado → SelectKBest no pipeline → RF; nested CV **agrupada por paciente** | Diagrama do pipeline + esquema dos folds |
| 11 | **Exploração: tecidos têm assinatura** | Normais agrupados por tecido; tumores dispersos | [tcga_pca.png] |
| 12 | **Resultados: tumor vs normal** | acc 0,984 · AUC 0,997; por tecido; estômago com 2 normais | [tcga_etapa1.png] + tabela por tecido |
| 13 | **Resultados: tecido de origem e hierárquico** | acc 0,980 (6 classes); hierárquico 0,971; confusão cólon ↔ estômago | [tcga_etapa2.png] / [tcga_hierarquico.png] |
| 14 | **O que o modelo aprendeu** | KANK1 (cólon), RBP5/LITAF (fígado), estômago intermediário | [tcga_heatmap_top_etapa2.png] |
| 15 | **"Bom demais?" Auditoria de vazamento** | CpGs aleatórios = mesmo resultado; rótulos embaralhados = 0,50; 1 CpG = 90%; por que um experimento "fácil" importa (CUP, controle positivo) | Tabela 4.5 |
| 16 | **Tentativa: plasma (GSE186575)** | Primeiro dataset escolhido (mais parecido com o artigo); 98 amostras, cfMeDIP; o EDA revela sexo, lote "mix" e idade | [eda_pca.png] + [eda_pca_autossomos.png] |
| 17 | **Por que o plasma falhou** | Detecção AUC 0,55; idade sozinha 0,70; 0 promotores com FDR; 526 por lote; demonstração de vazamento (0,84 com × 0,27 sem); decisão de troca | Tabela 5.4 + barras "câncer vs lote" |
| 18 | **Tecido × plasma × artigo** | Por que é fácil no tecido e difícil no sangue; o que o artigo faz para vencer isso | Tabela 6 |
| 19 | **Limitações e por que não é teste clínico** | Perguntas 14 e 16 resumidas | Lista |
| 20 | **Conclusão** | A assinatura existe na fonte; no plasma exige escala e tecnologia; prevenção de vazamento verificada; uso real (CUP) | Frase-síntese |

**Frase-síntese para o fechamento:** *"A assinatura de metilação dos tecidos é forte e fácil de capturar na fonte, o tumor, mas praticamente invisível no plasma com dados públicos e tecnologia de baixa resolução. Essa distância é exatamente o problema que o artigo resolve com escala e sequenciamento direcionado."*
