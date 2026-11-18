# Prompt para gerar os slides com IA

**Como usar:**
1. Anexe as 5 figuras de `results/figures/` listadas no fim do prompt.
2. Copie e cole o bloco "PROMPT" abaixo.
3. Depois de gerado, confira **todos os números** contra o plano (seção 4). A IA pode arredondar ou inventar valores.

---

## PROMPT

Crie uma apresentação acadêmica em **português do Brasil**, com **11 slides + 2 de apêndice**, para um seminário de mestrado da disciplina Genômica Computacional, com duração de 12–15 minutos. A plateia já domina metilação do DNA, CpG e cfDNA: **não crie slides explicando esses fundamentos**.

### Regras obrigatórias
- **Use exatamente os números que estão abaixo.** Não arredonde de outro jeito, não invente métricas e não crie gráficos com dados fictícios.
- **Gráficos de resultados:** use **somente as imagens anexadas**, com o nome de arquivo indicado em cada slide. Não recrie esses gráficos.
- **Diagramas e fluxogramas:** desenhe você mesmo quando eu pedir, com formas simples (caixas e setas), sem ícones decorativos em excesso.
- **Estilo:** limpo e científico. Fundo claro, uma cor de destaque (azul petróleo), fonte sem serifa e no máximo ~40 palavras de texto por slide. Os números-chave ficam grandes e em destaque.
- **Notas do apresentador:** em cada slide, 2–4 frases com o que falar.
- **Formato:** 16:9, exportável para PowerPoint e PDF.

---

### Slide 1 — Título
- **Título:** Detecção de câncer e tecido de origem por metilação do DNA
- **Subtítulo:** Um classificador com dados públicos, inspirado em Liu et al. (2020), *Annals of Oncology*
- **Rodapé:** João Pedro Thurow Tomaszewski · Genômica Computacional · 2026

### Slide 2 — O artigo: de onde vêm os dados
Layout em duas colunas.

**Coluna esquerda: estudos e amostras.**
- **CCGA** (NCT02889978): estudo caso-controle multicêntrico, com 142 centros. Forneceu todos os casos de câncer e parte dos controles.
- **STRIVE** (NCT03085888): mulheres em mamografia de rastreio. Forneceu controles extras, para medir especificidade acima de 99%.
- Tabela pequena:

| | Câncer | Sem câncer | Total |
|---|---|---|---|
| Treino | 1.742 | 2.978 | 4.720 |
| Validação | 740 | 1.229 | 1.969 |
| **Total** | **2.482** | **4.207** | **6.689** |

**Coluna direita: tecnologia.** Desenhe um fluxo horizontal de 4 caixas:
> Plasma (cfDNA) → Conversão por bissulfito → Captura direcionada (painel) → Sequenciamento (NovaSeq)

Embaixo do fluxo: "Lê o estado de cada CpG **em cada fragmento**".

Destaque: **103.456 regiões · ~1,1 milhão de CpGs**, escolhidas com WGBS do 1º sub-estudo + arrays de metilação do TCGA.

### Slide 3 — O artigo: classificação e vazamento
Desenhe um diagrama de **duas etapas**:
> Amostra → **(i) Detecção:** câncer vs não-câncer → se positivo → **(ii) Tecido de origem:** ~20 sítios

Embaixo do diagrama: "*Source models* por região → features → par de regressões logísticas".

Três números grandes:
- **99,3%** especificidade;
- **18% → 93%** sensibilidade do estádio I ao IV;
- **93%** acurácia do tecido de origem.

Caixa lateral "Contra vazamento":
- validação independente, cega e com o classificador travado antes;
- riscos: os controles do STRIVE são só mulheres de um centro; há várias amostras por pessoa.

### Slide 4 — Dados públicos: o que dá para usar?
- **Topo:** "CCGA/STRIVE: dados privados (GRAIL), sem acesso público (`ipdSharing: NO`)".
- Tabela compacta:

| Dataset | Material | Tecnologia | n | Situação |
|---|---|---|---|---|
| **TCGA (Xena)** | Tecido | Array 450K (bissulfito) | ~100–900/tipo | **Escolhido** |
| GSE186575 | Plasma | cfMeDIP-seq | 98 | Tentado |
| GSE81314 | Plasma | 5hmC (outra marca) | 64 | Descartado |
| GSE314062 | Plasma | Bissulfito direcionado | 5 | Pequeno demais |

- **Por que TCGA:** mesma química do artigo (bissulfito, por CpG), 6 tecidos, tecido normal pareado e ID do paciente (permite separar indivíduos).
- **Aviso em destaque (caixa amarela):** "Tecido ≠ plasma: o problema passa a ser **classificar o tecido**, não detectar câncer no sangue."

### Slide 5 — Tentativa: plasma (GSE186575)
- **Subtítulo:** "O dataset mais parecido com o artigo… e por que desistimos". 98 amostras de plasma: colorretal, gástrico e saudáveis.
- Três números grandes lado a lado, com legenda curta:
  - **AUC 0,55**: detecção de câncer no nível do acaso.
  - **0 × 526**: promotores significativos para câncer × para o lote técnico "mix".
  - **0,70**: AUC usando **só a idade** (melhor que a metilação).
- **Rodapé:** "Lote e idade dominam o sinal → trocamos para o TCGA. A tentativa está reproduzida no notebook."

### Slide 6 — TCGA: dados e pipeline sem vazamento
**Parte de cima: dados.**
- 1.290 amostras de 919 pacientes, 6 tecidos (cólon, estômago, pulmão, fígado, mama, rim).
- 380 normais e 910 tumores (~150 por tecido).
- 336 mil CpGs (sem valores ausentes, só autossomos) → 20 mil de maior variância.

**Parte de baixo: fluxograma do pipeline.**
> Matriz (amostras × CpGs) → [dentro de cada fold de treino: **SelectKBest (teste F)** → **Random Forest** (500 árvores, classes balanceadas)] → previsão no fold de teste

Desenhe uma moldura destacada em volta do bloco entre colchetes, com o rótulo "ajustado só no treino".

Ao lado, um esquema de folds: 5 blocos de pacientes, um deles marcado como teste. Legenda: "**Nested CV agrupada por paciente**: tumor e normal da mesma pessoa nunca ficam em treino e teste ao mesmo tempo (371 pacientes têm os dois)".

### Slide 7 — Resultados
- **Imagens:** `tcga_etapa1.png` (esquerda) e `tcga_etapa2.png` (direita).
- Números em destaque acima das imagens:
  - **Tumor vs normal:** acurácia **0,984** · ROC-AUC **0,997**
  - **Tecido de origem (6 classes):** acurácia **0,980** · F1 macro **0,980** · ROC-AUC **0,999**
  - **Fluxo hierárquico (normal + 6 tecidos):** acurácia **0,971**
- Linha de observação: "Pior caso: cólon (sensibilidade 0,94), confundido com estômago. Estômago tem só 2 normais."

### Slide 8 — "Bom demais?" Auditoria de vazamento
Tabela enxuta:

| Teste | Tumor vs normal | Tecido de origem |
|---|---|---|
| Referência | acc 0,982 | acc 0,978 |
| Pacientes em treino e teste | 0 | 0 |
| Filtro de variância dentro do fold | acc 0,982 | acc 0,977 |
| **Rótulos embaralhados** | **AUC 0,50** | **AUC 0,50** |
| **20 mil CpGs aleatórios** | acc 0,988 | acc 0,980 |
| 1 único CpG | acc 0,905 | acc 0,331 |

- **Caixa lateral "Contraprova no plasma":** com rótulos aleatórios, selecionar features **antes** da validação cruzada fabrica **AUC 0,84**; selecionando dentro dos folds, o resultado é **0,27**.
- **Conclusão em destaque:** "Sem vazamento: a tarefa em tecido é genuinamente fácil — o sinal está espalhado pelo genoma."

### Slide 9 — O que separa as classes
- **Imagem:** `tcga_heatmap_top_etapa2.png`, grande, à esquerda.
- À direita, 3 tópicos curtos:
  - **KANK1:** metilado só no cólon.
  - **RBP5 / LITAF:** desmetilados no fígado.
  - **Estômago intermediário** entre o cólon e os demais → explica a confusão cólon ↔ estômago.
- **Rodapé:** "Muitos CpGs redundantes: nenhum sozinho é indispensável."

### Slide 10 — Limites: por que isto não é um teste clínico
Três colunas com títulos:
- **Limitações:**
  - p ≫ n (20 mil features para 1.290 amostras);
  - estômago com 2 normais;
  - desbalanceamento tratado com pesos;
  - possível efeito de lote no TCGA.
- **Diferenças em relação ao artigo:**
  - tecido × cfDNA;
  - 1.290 × 6.689 amostras;
  - Random Forest × regressão logística;
  - validação cruzada × validação cega e travada;
  - limiar 0,5 × especificidade fixa de 99%.
- **Não é teste clínico:**
  - exige biópsia;
  - desenho caso-controle;
  - sem validação externa;
  - sem análise por estádio;
  - especificidade não medida na escala necessária.

### Slide 11 — Conclusão
Frase central grande:
> "A assinatura de metilação dos tecidos é forte e fácil de capturar na fonte, o tumor, mas praticamente invisível no plasma com dados públicos. Essa distância é exatamente o problema que o artigo resolve com escala e sequenciamento direcionado."

Embaixo, em uma linha: "Uso real em tecido: identificar a origem de câncer de sítio primário desconhecido (3–5% dos casos)."

### Apêndice A1 — PCA
- **Imagem:** `tcga_pca.png`.
- **Legenda:** "Normais agrupados por tecido; tumores dispersos."

### Apêndice A2 — Fluxo hierárquico
- **Imagem:** `tcga_hierarquico.png`.
- **Legenda:** "372 de 380 normais corretos; principal erro: 9 tumores de cólon previstos como estômago."

**Imagens anexadas:** `tcga_etapa1.png`, `tcga_etapa2.png`, `tcga_heatmap_top_etapa2.png`, `tcga_pca.png`, `tcga_hierarquico.png`.
