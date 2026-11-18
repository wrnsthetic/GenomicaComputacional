# Detecção de câncer e inferência do tecido de origem por padrões de metilação do DNA em cfDNA plasmático

**Disciplina:** Genômica Computacional
**Autor:** João Pedro Thurow Tomaszewski
**Artigo de referência:** Liu, M. C. et al. (2020). *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*. Annals of Oncology, 31(6), 745–759.
**Notebook reproduzível:** `notebooks/trabalho_final.ipynb`

> ⚠️ RASCUNHO — trechos marcados com `[[PREENCHER]]` dependem dos resultados finais do notebook.

---

## Resumo

Este trabalho estuda a metodologia de detecção multicâncer baseada em metilação de DNA livre circulante (cfDNA) proposta por Liu et al. (2020) e implementa uma versão experimental, em pequena escala, com dados públicos. Como os dados dos estudos CCGA e STRIVE não estão publicamente disponíveis, investigamos conjuntos alternativos no NCBI GEO e no TCGA e selecionamos o GSE186575, composto por perfis de cfMeDIP-seq do plasma de 32 pacientes com câncer colorretal, 35 com câncer gástrico e 31 indivíduos saudáveis. Replicamos a estrutura em duas etapas do artigo, com detecção de câncer seguida da predição do tecido de origem, usando Random Forest e validação cruzada aninhada, com todas as etapas de seleção de features restritas ao treino. A análise exploratória revelou dois fatores de confundimento relevantes: o sexo biológico, removido ao excluir os cromossomos sexuais, e um provável efeito de lote associado aos controles. Esses fatores são discutidos como limitações. `[[PREENCHER: principais números de desempenho]]`

---

## 1. Introdução

O câncer é diagnosticado com frequência em estágios avançados, quando as opções terapêuticas são mais limitadas e a mortalidade é maior. Programas de rastreamento eficazes existem para poucos tipos de tumor (mama, colo do útero, colorretal e pulmão em populações de risco), e cada um deles mira um único órgão. Um exame de sangue capaz de detectar vários tipos de câncer ao mesmo tempo e de indicar o órgão de origem mudaria esse cenário. É essa a proposta do teste avaliado no artigo de referência.

Os objetivos deste trabalho são:
1. compreender criticamente a metodologia do artigo (Parte I);
2. investigar a disponibilidade de dados públicos de metilação adequados ao problema (Parte II);
3. implementar e avaliar um classificador experimental de metilação com dados acessíveis (Parte III);
4. discutir os resultados, as limitações e as diferenças em relação ao trabalho original (Parte IV).

---

## 2. Parte I — Estudo do artigo

### 2.1 Fundamentação biológica

**Metilação do DNA.** É uma modificação epigenética: altera a regulação dos genes sem alterar a sequência de bases. Em mamíferos, consiste principalmente na adição de um grupo metil ao carbono 5 da citosina (5-metilcitosina, 5mC), quase sempre em dinucleotídeos CpG. O genoma humano tem cerca de 28 a 30 milhões de CpGs. A maioria está metilada, exceto nas *ilhas CpG*, regiões ricas em CpG geralmente localizadas em promotores de genes. A metilação de um promotor tende a silenciar o gene correspondente. O padrão de metilação é estabelecido durante a diferenciação celular e é estável e específico de cada tipo celular. Por isso funciona como uma "impressão digital" do tecido.

**Metilação e câncer.** Os tumores apresentam duas alterações características de metilação:
- *hipometilação global*, associada a instabilidade genômica e à ativação de elementos repetitivos;
- *hipermetilação focal* de ilhas CpG em promotores, que pode silenciar genes supressores de tumor e genes de reparo de DNA.

Essas alterações surgem cedo na carcinogênese e se repetem entre pacientes, o que as torna bons biomarcadores.

**DNA livre circulante (cfDNA) e biópsia líquida.** O cfDNA é formado por fragmentos curtos de DNA (cerca de 167 pb, o tamanho de um nucleossomo) liberados no sangue principalmente por morte celular (apoptose e necrose). Em pessoas saudáveis, a maior parte vem de células hematopoéticas: leucócitos e precursores eritroides (Moss et al., 2018). Em pacientes com câncer, uma fração, o *DNA tumoral circulante* (ctDNA), vem do tumor. Essa fração costuma ser baixa (muitas vezes menor que 1%), principalmente em estágios iniciais. A *biópsia líquida* é a análise desse material a partir de uma simples coleta de sangue, sem procedimento invasivo no órgão.

**Por que a metilação permite detectar o câncer e inferir sua origem.** O ctDNA carrega o padrão de metilação da célula que o liberou. Isso permite responder a duas perguntas:
1. *Existe DNA com padrão tumoral no plasma?* Fragmentos com padrões aberrantes, ausentes no cfDNA de pessoas saudáveis, indicam a presença de câncer (detecção).
2. *De qual tecido ele vem?* Como o padrão de metilação também reflete o tecido de origem da célula tumoral, é possível inferir o órgão afetado (*tissue of origin*, TOO).

Mutações, por outro lado, raramente indicam o tecido de origem. Além disso, são confundidas pela hematopoiese clonal de potencial indeterminado (CHIP), que gera mutações em leucócitos de pessoas sem câncer.

### 2.2 Caracterização das amostras e dos dados

**Os estudos CCGA e STRIVE.**
- **CCGA** (*Circulating Cell-free Genome Atlas*, NCT02889978) é um estudo prospectivo, multicêntrico, observacional e do tipo caso-controle, com acompanhamento longitudinal. Coletou amostras de 15.254 participantes (56% com câncer e 44% sem câncer) em 142 centros da América do Norte. Foi dividido em três sub-estudos pré-especificados: (1) descoberta, que comparou WGBS, WGS e sequenciamento direcionado de mutações e elegeu a metilação como a melhor abordagem; (2) treino e validação do ensaio de metilação direcionada, que é o artigo estudado; (3) validação adicional. O CCGA forneceu todos os casos de câncer e parte dos controles.
- **STRIVE** (NCT03085888) é uma coorte prospectiva de 99.286 mulheres submetidas a mamografia de rastreamento em 35 centros. Neste trabalho, o STRIVE forneceu **controles adicionais sem câncer**, todos de um único centro. Eles eram necessários porque, para operar com especificidade acima de 99%, é preciso um grande número de controles, tanto para treinar o classificador quanto para estimar a especificidade com precisão (o objetivo era especificidade acima de 99% com mais de 90% de confiança).

**Composição dos conjuntos.**

| Conjunto | Com câncer | Sem câncer (CCGA + STRIVE) | Total |
|---|---|---|---|
| Treinamento | 1.742 | 2.978 (1.391 + 1.587) | 4.720 |
| Validação | 740 | 1.229 (614 + 615) | 1.969 |
| **Total** | **2.482** | **4.207** | **6.689** |

Após critérios de elegibilidade e avaliabilidade, 4.316 participantes compuseram a população da análise primária (3.052 no treino e 1.264 na validação). Outras 354 amostras do CCGA foram reservadas para um conjunto de referência de biópsias tumorais. A distribuição de estádios foi semelhante entre treino e validação (cerca de 28% no estádio I, 25% no II, 20–21% no III e 23–24% no IV).

**Tecnologia.** O cfDNA plasmático (até 75 ng) foi submetido a conversão por bissulfito, que transforma citosinas não metiladas em uracila e preserva as metiladas. Em seguida, foram feitos o preparo de bibliotecas com índice duplo, o enriquecimento por captura híbrida com um painel customizado (Twist Bioscience) e o sequenciamento *paired-end* de 150 pb em Illumina NovaSeq. A profundidade mediana foi de 113 milhões de fragmentos (139× de cobertura única no alvo). A leitura por bissulfito permite observar o estado de metilação de cada CpG **em cada fragmento**, ou seja, o padrão de metilação molécula a molécula.

**Seleção das regiões genômicas.** As regiões do painel foram escolhidas a partir de duas fontes:
- dados de WGBS do primeiro sub-estudo (cfDNA de pacientes e controles, biópsias tumorais, leucócitos e tecidos comerciais);
- arrays de metilação do TCGA.

A partir delas, foram identificadas regiões do genoma (hg19) com padrões de metilação específicos de câncer ou de tecido em comparação com o cfDNA de controles. Um algoritmo customizado combinou os alvos mais informativos. O painel final cobre **103.456 regiões (17,2 Mb) e cerca de 1,1 milhão de CpGs**, uma fração pequena e informativa dos cerca de 30 milhões de CpGs do genoma. Isso permite sequenciar essas regiões com muito mais profundidade do que no WGBS, pelo mesmo custo.

**Separação dos dados, validação independente e riscos de vazamento.** O desenho do estudo tem várias salvaguardas:
- **Conjunto de validação independente**, mantido separado durante todo o desenvolvimento;
- **Classificador travado** (*locked*) e escolhido antes da liberação dos dados de validação;
- **Análise duplo-cega:** os desenvolvedores não conheciam o status de câncer da validação, e uma equipe independente de integridade de dados revisou a base antes da liberação;
- **Validação cruzada no treino** para ajustar os dois níveis de modelo e os limiares, garantindo que os dados dos *folds* de teste não fossem acessados durante o treino;
- **Distribuição pré-especificada** de tipos de câncer, estádios e centros entre os sub-estudos, com dados demográficos comparáveis entre casos e controles;
- **Análise de sensibilidade por centro** (*site balancing*): o desempenho em centros removidos do treino foi comparável (sensibilidade de 53,6% contra 50,0%);
- **Comparação da especificidade por origem dos controles:** taxa de falsos positivos de 0,7% (CCGA) contra 0,6% (STRIVE), sugerindo que o classificador não aprendeu a origem do controle.

Ainda assim, há riscos residuais que merecem atenção:
1. Os controles do STRIVE são **todas mulheres de um único centro**. Isso cria a possibilidade de confundimento por sexo, centro ou lote de processamento. Os autores avaliam esse ponto pela comparação de taxas de falsos positivos.
2. As regiões do painel foram escolhidas com dados do primeiro sub-estudo. Se houvesse sobreposição de participantes com a validação, isso seria vazamento. O artigo informa que as biópsias do segundo sub-estudo usadas no WGBS não participaram da avaliação do classificador.
3. Os modelos foram treinados em 11.154 amostras de 5.854 participantes, ou seja, **mais de uma amostra por pessoa**. Sem divisão por indivíduo, amostras da mesma pessoa poderiam cair em treino e teste.
4. Nem todos os controles tinham 1 ano completo de acompanhamento, então alguns "não-câncer" podem ter câncer ainda não diagnosticado.

### 2.3 Estratégia de classificação

O artigo resolve **dois problemas encadeados**:

| | (i) Detecção | (ii) Tecido de origem (TOO) |
|---|---|---|
| Pergunta | Há sinal de câncer no plasma? | De qual órgão vem o sinal? |
| Amostras | Todas | Apenas as positivas em (i) |
| Saída | Câncer / não-câncer | Um entre cerca de 20 sítios/tipos tumorais |
| Métrica central | Especificidade (fixada em ~99%) e sensibilidade | Acurácia da predição de TOO |

**Arquitetura.** Primeiro, *source models* analisam os fragmentos de cada região e calculam o quanto o padrão de metilação observado se parece com o de cada tipo de câncer, em comparação com o cfDNA não canceroso. As saídas desses modelos formam as features da etapa seguinte. Depois, um **par de regressões logísticas em *ensemble*** faz a classificação: uma decide câncer ou não-câncer, e a outra atribui o TOO. O limiar de decisão foi ajustado para especificidade de cerca de 99%, porque, num rastreamento populacional, falsos positivos geram exames invasivos desnecessários.

**Resultados principais (validação).**
- Especificidade de 99,3% (IC 95%: 98,3–99,8%), ou seja, taxa de falsos positivos de 0,7%.
- Sensibilidade nos estádios I a III: 67,3% nos 12 cânceres pré-especificados e 43,9% em todos os tipos.
- A sensibilidade cresce com o estádio: de 18% (I) para 43% (II), 81% (III) e 93% (IV), considerando todos os tipos.
- O TOO foi previsto em 96% das amostras com sinal de câncer, com 93% de acurácia.

---

## 3. Parte II — Busca e avaliação de datasets públicos

### 3.1 Disponibilidade dos dados CCGA e STRIVE

- O artigo **não contém declaração de disponibilidade de dados** nem número de acesso em repositório (GEO, EGA, dbGaP).
- No ClinicalTrials.gov, o registro do **STRIVE (NCT03085888) declara `IPD sharing: No`**, isto é, não há compartilhamento de dados individuais. O registro do **CCGA (NCT02889978)** não declara nenhum plano de compartilhamento.
- O estudo foi financiado pela GRAIL, Inc., que comercializa o teste (Galleri). Os dados são, portanto, **proprietários**.

Para ilustrar a diferença entre dados públicos e dados de **acesso controlado**: o maior compêndio de cfMeDIP-seq pan-câncer recente (Zeng et al., 2026; 1.294 amostras de 11 tipos de câncer) deposita seus dados no EGA. O acesso exige um acordo formal (*Data Access Agreement*) avaliado por um comitê institucional, um processo incompatível com o prazo da disciplina.

### 3.2 Conjuntos candidatos

Todos os metadados foram conferidos diretamente nos registros do GEO e do UCSC Xena em 17/11/2026.

| Campo | GSE186575 | GSE81314 | TCGA (Xena) | GSE314062 |
|---|---|---|---|---|
| Publicação | PMID 37992683 | Li et al., 2017 (*Cell Res.*) | TCGA / PanCanAtlas | PMID 42524278 |
| Nº de amostras | 98 | 126 GSMs; matriz com 64 pacientes | ~100–900 por tipo tumoral | 5 |
| Classes | CRC 32, gástrico 35, saudável 31 | Pulmão 15, fígado 10, pâncreas 7, estômago 5, cólon 4, mama 4, GBM 4 | Tumor vs normal adjacente, por tipo | 3 pacientes HNC |
| Controles | 31 saudáveis | 8 saudáveis + 7 HBV | Tecido normal adjacente | Nenhum |
| Material | **Plasma (cfDNA)** | Plasma (cfDNA) | Tecido tumoral | Plasma |
| Tecnologia | cfMeDIP-seq | hMe-Seal (5hmC) | Illumina 450K | Bissulfito direcionado |
| Formato | Matriz CPM por promotor (58.288 × 98) | Matriz por gene (20.244 × 68) | Beta values por CpG (~485 mil) | Arquivos brutos |
| Acesso | Público | Público | Público | Público |
| Viabilidade | **Alta** | Média | Alta (só tecido) | Inviável |

Também foram avaliados o GSE275598 (cfMeDIP-seq de tumores cerebrais pediátricos, fora do escopo multicâncer em adultos) e o GSE122126 (atlas de metilação de tipos celulares, feito para deconvolução e com pouquíssimos cânceres).

**Observações sobre cada candidato:**
- **GSE81314:** mede 5-hidroximetilcitosina (5hmC), uma marca epigenética **diferente** da 5mC usada no artigo. Tem de 4 a 15 pacientes por classe, pouco para um classificador multiclasse. A matriz também contém colunas que não são pacientes (DNA genômico de leucócitos e controles técnicos de *input*), que precisariam ser removidas.
- **TCGA:** grande e bem anotado, mas formado por **tecido tumoral**, não por plasma. Classificar tecido tumoral é um problema mais fácil e diferente da detecção em cfDNA, onde o sinal tumoral é diluído em DNA de células sanguíneas. O arquivo pan-câncer único tem 14 GB, inviável no Google Colab gratuito. Os arquivos por coorte (100–800 MB) são utilizáveis.

### 3.3 Dataset escolhido e justificativa

Escolhemos o **GSE186575** por ser o único candidato público que reúne ao mesmo tempo:
1. **cfDNA plasmático**, o mesmo material biológico do artigo;
2. **controles saudáveis**, necessários para a tarefa de detecção;
3. **mais de um tipo de câncer**, necessário para a tarefa de tecido de origem;
4. tamanho amostral razoável por classe (31 a 35);
5. dados já processados em matriz, leves (22 MB) e de acesso livre.

Assim, é possível reproduzir **em pequena escala o fluxo em duas etapas do artigo**: detecção seguida de TOO. O TCGA fica registrado como alternativa para um experimento complementar de classificação de tecido.

---

## 4. Parte III — Metodologia

### 4.1 Dados

O **GSE186575** reúne amostras de plasma coletadas no Peking University First Hospital (China). O cfDNA foi extraído de 1 a 3 mL de plasma e perfilado por **cfMeDIP-seq** (*cell-free methylated DNA immunoprecipitation*). Nessa técnica, um anticorpo contra 5mC "pesca" os fragmentos metilados, que são então sequenciados. Os autores alinharam as leituras ao genoma hg38 e contaram os fragmentos por **região promotora** de cada gene, em contagens por milhão (CPM).

- **Matriz:** 98 amostras × 58.288 promotores (CPM). Não há valores ausentes; cerca de 9,3% dos valores são zero.
- **Rótulos:** derivados do prefixo de cada amostra: `CRC` (colorretal, n = 32), `STAD` (adenocarcinoma gástrico, n = 35) e `NC` (controle saudável, n = 31).
- **Metadados:** sexo, idade e estadiamento, extraídos da *series matrix* do GEO.
- **Indivíduos:** cada amostra corresponde a um participante distinto, então não há múltiplas amostras por pessoa.

**Diferença importante em relação ao artigo:** no cfMeDIP-seq, a feature é a **quantidade de DNA metilado em uma região** (sinal agregado por promotor). No artigo, a medida é o **padrão de metilação de cada CpG em cada fragmento**, obtido por bissulfito.

### 4.2 Análise exploratória e pré-processamento

A análise exploratória por PCA, feita nos 5.000 promotores de maior variância, revelou três fatores não biológicos com impacto potencial no classificador.

1. **Sexo biológico.** Com todos os cromossomos incluídos, o 2º componente principal (12,7% da variância) separava as amostras quase perfeitamente por sexo. Isso se explica pelos promotores dos cromossomos X e Y. O fator é crítico porque o grupo gástrico é 77% masculino, enquanto os demais são balanceados: um classificador poderia distinguir STAD de CRC pelo sexo. **Decisão:** remover os promotores dos cromossomos X, Y, mitocondrial e os sem anotação (mapeamento gene → cromossomo via Ensembl BioMart, hg38), restando **53.768 promotores autossômicos**. Depois da remoção, o componente associado ao sexo desaparece.

2. **Componente técnico dominante.** Mesmo só com os autossomos, o PC1 explica 46% da variância e está correlacionado com a fração de zeros por amostra (r = −0,63), não com o grupo. Isso sugere variação na eficiência da imunoprecipitação ou na profundidade de sequenciamento. Duas amostras de CRC aparecem como outliers extremos nesse eixo.

3. **Provável efeito de lote.** 27 amostras têm o sufixo `mix` no nome, e elas se distribuem de forma muito desigual: 25 dos 31 controles, 2 dos 32 CRC e nenhum dos 35 STAD. O significado do sufixo não é documentado no GEO, então o tratamos como **hipótese** de lote de processamento. Um Random Forest consegue distinguir controles `mix` de controles não-`mix` (todos saudáveis) com AUC de cerca de 0,82, o que confirma uma assinatura técnica associada a esse grupo.

Além disso, os **controles são cerca de 10 anos mais jovens** (mediana de 50 anos contra 62–64 nos grupos com câncer). A idade é um conhecido modulador da metilação do DNA.

**Pré-processamento final:**
- transformação `log2(CPM + 1)`, aplicada amostra a amostra, para reduzir a assimetria;
- remoção dos cromossomos sexuais e mitocondrial.

Nenhuma dessas etapas usa estatísticas de outras amostras, então não há risco de vazamento.

### 4.3 Tratamento dos fatores de confundimento

Sexo, lote e idade estão correlacionados com as classes. Por isso, o classificador pode aprender esses fatores em vez do sinal de metilação associado ao câncer. Não é possível eliminá-los totalmente com os dados disponíveis, então as estratégias foram:
- **Sexo:** remoção das features que carregam o sinal (promotores de chrX/Y).
- **Lote `mix`:** análise de sensibilidade nas amostras que "cruzam" o lote, ou seja, os 6 controles não-`mix` e os 2 CRC `mix`. Se o modelo tivesse aprendido o lote, erraria justamente nessas amostras.
- **Idade:** um modelo usando apenas a idade serve de linha de base.

Métodos de correção de lote, como o ComBat (Johnson et al., 2007), não foram aplicados. Com o lote quase alinhado à classe, a correção removeria também o sinal biológico, já que não há como separar o efeito do lote do efeito do câncer (Leek et al., 2010).

### 4.4 Classificador e validação

**Algoritmo: Random Forest** (Breiman, 2001). É um conjunto de árvores de decisão, cada uma treinada numa amostra *bootstrap* dos dados e considerando, a cada divisão, um subconjunto aleatório de features. A predição final é a votação das árvores. O Random Forest foi escolhido porque:
- lida bem com alta dimensionalidade (p ≫ n);
- captura interações e relações não lineares;
- é pouco sensível à escala das features;
- fornece medidas de importância das features.

Usamos 500 árvores e `class_weight="balanced"` para compensar o desbalanceamento entre classes.

**Pipeline** (todas as etapas ajustadas **somente no treino** de cada partição):
1. **Filtro de presença:** mantém promotores com sinal > 0 em pelo menos 50% das amostras de treino.
2. **Seleção supervisionada:** `SelectKBest` com teste F (ANOVA), mantendo os *k* promotores mais associados ao rótulo.
3. **Random Forest.**

**Validação:**
- **Validação cruzada aninhada** (*nested CV*; Varma & Simon, 2006). O laço externo usa 5 *folds* estratificados repetidos 3 vezes e estima o desempenho. O laço interno usa 5 *folds* e escolhe os hiperparâmetros: *k* ∈ {100, 500, 2000} e `max_features` ∈ {√p, 0,1·p}.
- **Hold-out estratificado (70/30)** para apresentar a matriz de confusão, a curva ROC e a importância das features num conjunto nunca usado no treino nem na escolha de hiperparâmetros.

**Tarefas:**
1. **Detecção:** câncer (CRC + STAD, n = 67) vs saudável (n = 31).
2. **TOO:** CRC (n = 32) vs STAD (n = 35), só com amostras de câncer.
3. **Fluxo hierárquico:** predições fora da amostra (*out-of-fold*) da etapa 1 e, nos positivos, da etapa 2, resultando numa classificação final em três classes (NC / CRC / STAD).

**Métricas:** matriz de confusão, precisão, sensibilidade (*recall*), especificidade, F1-score, acurácia e ROC-AUC.

**Importância das features:** *permutation importance* no conjunto de teste, que mede a queda de AUC ao embaralhar cada promotor. Os promotores são anotados com o nome do gene via Ensembl.

---

## 5. Resultados

`[[PREENCHER com as saídas do notebook]]`

### 5.1 Detecção: câncer vs saudável
### 5.2 Tecido de origem: CRC vs STAD
### 5.3 Fluxo hierárquico
### 5.4 Confundidores: lote e idade
### 5.5 Demonstração do efeito do vazamento
### 5.6 Promotores mais relevantes

---

## 6. Parte IV — Discussão

**6.1 Há diferenças de metilação capazes de separar as classes?** `[[PREENCHER]]`

**6.2 Quais regiões contribuem mais?** `[[PREENCHER: genes do top de importância e relação com a literatura]]`

**6.3 Desempenho e variação entre classes.** `[[PREENCHER]]`

**6.4 Limitações: tamanho amostral, desbalanceamento e alta dimensionalidade.**
- Com 98 amostras e mais de 53 mil features (p ≫ n), há alto risco de sobreajuste e de associações espúrias. Isso foi mitigado pela seleção de features dentro da validação e pela CV aninhada, mas não eliminado.
- Os intervalos de confiança são largos: com 31 controles, uma especificidade estimada em 100% é compatível, com 95% de confiança, com valores reais em torno de 89%. É impossível verificar a especificidade de 99% exigida para rastreamento.
- O desbalanceamento é moderado (cerca de 2:1 na detecção), tratado com `class_weight="balanced"` e métricas além da acurácia.
- Os dados vêm de um único centro, com possível lote confundido com a classe e controles mais jovens.

**6.5 Diferenças em relação à metodologia original.**

| Aspecto | Liu et al. (2020) | Este trabalho |
|---|---|---|
| Amostras | 6.689 (multicêntrico, EUA) | 98 (um centro, China) |
| Tipos de câncer | > 50 | 2 (colorretal e gástrico) |
| Tecnologia | Bissulfito + captura direcionada (1,1 mi CpGs) | cfMeDIP-seq (imunoprecipitação) |
| Resolução | CpG por fragmento | Sinal agregado por promotor |
| Features | *Source models* por região | Promotores selecionados por teste F |
| Classificador | Regressões logísticas em *ensemble* | Random Forest |
| Limiar | Ajustado para especificidade de ~99% | 0,5 (padrão) |
| Validação | Conjunto independente, cego, classificador travado | CV aninhada + hold-out |
| Controles de viés | Balanceamento de centros, comparação CCGA × STRIVE | Remoção de chrX/Y, análise de lote e idade |

**6.6 Por que o resultado não pode ser interpretado como um teste clínico de detecção precoce.**
1. **Desenho caso-controle:** os pacientes já tinham diagnóstico e muitos eram sintomáticos. Um teste de rastreamento precisa ser avaliado numa população assintomática, onde a prevalência é baixa (cerca de 1,3% ao ano) e o valor preditivo positivo depende fortemente da especificidade.
2. **Especificidade não mensurável:** com 31 controles, não é possível demonstrar uma taxa de falsos positivos abaixo de 1%.
3. **Ausência de validação externa:** não houve conjunto independente de outro centro ou lote, e há evidência de confundimento por lote e idade.
4. **Só dois tipos de câncer**, ambos gastrointestinais. O teste não foi avaliado contra outros cânceres nem contra doenças benignas que alteram o cfDNA (inflamação, por exemplo).
5. **Sem estratificação por estádio** com poder suficiente: a detecção precoce depende do desempenho no estádio I, onde a fração tumoral é mínima.
6. **Sem limiar pré-especificado nem classificador travado** antes da avaliação.
7. Não há demonstração de **benefício clínico**, como redução de mortalidade, nem análise de custos e danos.

---

## 7. Conclusão

`[[PREENCHER]]`

---

## Registro das fontes, classes e validação

| Item | Descrição |
|---|---|
| Dados de metilação | GEO GSE186575: `GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz` e `GSE186575_series_matrix.txt.gz` (FTP público do NCBI) |
| Anotação de genes | Ensembl BioMart (GRCh38): ID do gene, nome, cromossomo, biotipo |
| Classes | NC = saudável (31); CRC = câncer colorretal (32); STAD = câncer gástrico (35). Detecção: CRC+STAD vs NC. TOO: CRC vs STAD |
| Validação | Nested CV (externa: 5 folds estratificados × 3 repetições; interna: 5 folds) + hold-out estratificado 70/30; seleção de features dentro do pipeline; uma amostra por indivíduo |

---

## Referências

1. Liu, M. C. et al. (2020). Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA. *Annals of Oncology*, 31(6), 745–759. doi:10.1016/j.annonc.2020.02.011
2. ClinicalTrials.gov. NCT02889978 (CCGA) e NCT03085888 (STRIVE).
3. NCBI GEO. GSE186575 — Plasma Multi-omics Nucleic Acid Landscape of Gastrointestinal Cancer [MeDIP-seq]. PMID 37992683.
4. Li, W. et al. (2017). 5-Hydroxymethylcytosine signatures in circulating cell-free DNA as diagnostic biomarkers for human cancers. *Cell Research*, 27, 1243–1257. (GSE81314)
5. Moss, J. et al. (2018). Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. *Nature Communications*, 9, 5068. (GSE122126)
6. Shen, S. Y. et al. (2018). Sensitive tumour detection and classification using plasma cell-free DNA methylomes. *Nature*, 563, 579–583.
7. Zeng et al. (2026). A pan-cancer compendium of 1,294 plasma cell-free DNA methylomes and fragmentomes enabling multicancer detection. *Nature Cancer*. doi:10.1038/s43018-026-01116-3
8. Goldman, M. J. et al. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology*, 38, 675–678.
9. Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32.
10. Varma, S. & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, 7, 91.
11. Ambroise, C. & McLachlan, G. J. (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS*, 99(10), 6562–6566.
12. Leek, J. T. et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*, 11, 733–739.
13. Johnson, W. E., Li, C. & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118–127.
14. Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825–2830.
