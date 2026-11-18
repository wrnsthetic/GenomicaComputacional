# Detecção de câncer e inferência do tecido de origem por padrões de metilação do DNA

**Disciplina:** Genômica Computacional
**Autor:** Laís Goulart
**Artigo de referência:** Liu, M. C. et al. (2020). *Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA*. Annals of Oncology, 31(6), 745–759.
**Notebook reproduzível:** `notebooks/trabalho_final.ipynb`

---

## Resumo

Este trabalho estuda a metodologia de detecção multicâncer baseada em metilação de DNA livre circulante (cfDNA) proposta por Liu et al. (2020) e implementa uma versão experimental, em pequena escala, com dados públicos. Os dados dos estudos CCGA e STRIVE não estão disponíveis publicamente, então avaliamos seis conjuntos alternativos.

A primeira tentativa usou o GSE186575, com perfis de cfMeDIP-seq do plasma de pacientes com câncer colorretal, câncer gástrico e indivíduos saudáveis. A detecção de câncer ficou no nível do acaso (ROC-AUC de 0,55). Nenhum promotor diferiu significativamente entre as classes, e o sinal era dominado por um efeito de lote e pela idade. Por isso, o experimento principal foi feito com dados de **tecido** do TCGA: perfis Illumina 450K de 1.290 amostras de 919 pacientes, cobrindo seis tipos de câncer (cólon, estômago, pulmão, fígado, mama e rim).

Replicamos a estrutura em duas etapas do artigo: primeiro tumor vs normal, depois o tecido de origem. O classificador é um Random Forest com seleção de features dentro do pipeline, avaliado por validação cruzada aninhada **agrupada por paciente**. Os resultados foram:
- tumor vs normal: acurácia de 0,984 e ROC-AUC de 0,997;
- tecido de origem (6 classes): acurácia de 0,980 e F1 macro de 0,980;
- fluxo hierárquico completo: acurácia de 0,971.

Uma auditoria descartou vazamento de informação. Com rótulos embaralhados, a AUC cai para 0,50, e 20 mil CpGs escolhidos ao acaso reproduzem o desempenho. A principal confusão foi entre cólon e estômago. Os CpGs mais informativos são específicos de tecido (KANK1 no cólon; RBP5 e LITAF no fígado).

O contraste entre tecido e plasma mostra que a assinatura de metilação é forte na fonte, mas fica diluída no sangue. Essa distância é justamente o que o artigo supera com escala e sequenciamento direcionado.

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
| Material | Plasma (cfDNA) | Plasma (cfDNA) | **Tecido** | Plasma |
| Tecnologia | cfMeDIP-seq | hMe-Seal (5hmC) | Illumina 450K (bissulfito) | Bissulfito direcionado |
| Formato | Matriz CPM por promotor (58.288 × 98) | Matriz por gene (20.244 × 68) | Beta values por CpG (~485 mil) | Arquivos brutos |
| Acesso | Público | Público | Público | Público |
| Viabilidade | Alta (**tentado**, seção 3.3) | Média | Alta (**escolhido**) | Inviável |

Também avaliamos dois outros conjuntos, ambos descartados:
- **GSE275598:** cfMeDIP-seq de tumores cerebrais pediátricos, fora do escopo multicâncer em adultos.
- **GSE122126:** atlas de metilação de tipos celulares, feito para deconvolução e com pouquíssimos cânceres.

**Observações sobre cada candidato:**
- **GSE81314:** mede 5-hidroximetilcitosina (5hmC), uma marca epigenética **diferente** da 5mC usada no artigo. Tem de 4 a 15 pacientes por classe, pouco para um classificador multiclasse.
- **GSE186575:** é o único candidato público que reúne cfDNA plasmático (o mesmo material do artigo), controles saudáveis e mais de um tipo de câncer. Por isso foi a primeira escolha.
- **TCGA:** grande, bem anotado, com vários tipos de câncer e com tecido normal pareado de muitos pacientes. A ressalva é que é **tecido tumoral**, não plasma. O arquivo pan-câncer único tem 14 GB, inviável no Google Colab gratuito, mas os arquivos por coorte (300–800 MB) são utilizáveis.

### 3.3 Avaliação experimental do GSE186575 (plasma)

Antes de adotar o GSE186575, avaliamos experimentalmente se ele sustentava um classificador. Esta seção resume essa avaliação, que motivou a troca de dataset.

**Dados.** O GSE186575 tem plasma do Peking University First Hospital (China), perfilado por **cfMeDIP-seq**. Nessa técnica, um anticorpo contra 5mC "pesca" os fragmentos metilados, que são então sequenciados e contados por **região promotora** de cada gene, em contagens por milhão (CPM). A matriz tem 98 amostras × 58.288 promotores, uma amostra por participante. A feature é a **quantidade de DNA metilado em uma região**, um sinal agregado. No artigo, a medida é o padrão de metilação de cada CpG em cada fragmento.

**Fatores de confundimento identificados na análise exploratória:**
1. **Sexo.** O 2º componente principal separava as amostras por sexo, por causa dos promotores dos cromossomos X e Y. O grupo gástrico é 77% masculino, então o classificador poderia distinguir os dois cânceres pelo sexo. Removemos os promotores de chrX, chrY, do DNA mitocondrial e os sem anotação (Ensembl BioMart, hg38). Restaram 53.768 promotores autossômicos.
2. **Provável efeito de lote.** 27 amostras têm o sufixo `mix` no nome, distribuídas de forma muito desigual: 25 dos 31 controles, 2 dos 32 CRC e nenhum dos 35 gástricos. O significado do sufixo não é documentado no GEO, então o tratamos como hipótese de lote de processamento.
3. **Idade.** Os controles são cerca de 10 anos mais jovens (mediana de 50 anos, contra 62–64 nos grupos com câncer).

**Modelagem.** Usamos a mesma lógica do experimento principal (seção 4.3):
- pipeline com filtro de presença (sinal em ≥ 50% das amostras de treino), `SelectKBest(f_classif)` e Random Forest;
- transformação `log2(CPM+1)`;
- validação cruzada aninhada, com laço externo de 5 *folds* × 3 repetições e laço interno de 5 *folds*, escolhendo *k* e `max_features`.

**Resultados.**

| Tarefa | ROC-AUC | Acurácia | Sensibilidade | Especificidade | F1 |
|---|---|---|---|---|---|
| Detecção: câncer (67) vs saudável (31) | 0,553 ± 0,182 | 0,623 | 0,763 | **0,324** | 0,728 |
| Origem: colorretal (32) vs gástrico (35) | 0,789 ± 0,136 | 0,711 | 0,752 | 0,663 | 0,726 |

A investigação da falha na detecção mostrou:
- **Nenhuma diferença significativa entre as classes.** No teste F promotor a promotor, com correção de Benjamini–Hochberg (FDR de 5%), houve **0 promotores significativos** em câncer vs saudável e 0 em colorretal vs gástrico. Já entre controles `mix` e não-`mix`, todos saudáveis, foram **526**. O lote deixa mais marca na metilação do que o câncer.
- **Lote:** os 6 controles não-`mix` foram **todos** classificados como câncer.
- **Idade:** um modelo usando apenas a idade atinge AUC de **0,70** na detecção, mais que o modelo de metilação.

Não aplicamos correção de lote (ComBat; Johnson et al., 2007). Com o lote quase alinhado à classe, a correção removeria também o sinal biológico (Leek et al., 2010).

**Demonstração do efeito do vazamento.** Com **rótulos aleatórios** (sem sinal real), selecionar os 500 promotores mais associados ao rótulo usando todas as amostras, antes da validação cruzada, produz AUC de **0,84**: um sinal inteiramente fabricado. Com a seleção dentro do pipeline, a AUC é de 0,27. O valor abaixo de 0,5 é flutuação esperada com n pequeno e rótulos sem sinal. Toda seleção supervisionada deste trabalho é feita dentro do pipeline.

**Conclusão da avaliação.** Detecção no nível do acaso, ausência de promotores significativos e lote mais forte que o câncer não permitem um experimento de classificação viável. A tentativa é mantida como evidência para a discussão: com dados públicos de pequena escala e baixa resolução, o sinal tumoral no plasma é superado por ruído técnico e demográfico.

### 3.4 Dataset escolhido e justificativa

O experimento principal usa o **TCGA**, com metilação de tecido (Illumina 450K) de seis tipos de câncer. Os motivos são:
1. **Mesma química do artigo:** conversão por bissulfito com leitura por CpG (*beta value*). O GSE186575 usa imunoprecipitação.
2. **Escala:** centenas de amostras por tipo de câncer, o que permite uma tarefa de tecido de origem com várias classes.
3. **Tecido normal adjacente** pareado de muitos pacientes, o que permite a etapa tumor vs normal.
4. **Identificação do paciente** no código da amostra, o que permite validação separando indivíduos.
5. **Relação direta com o artigo:** os autores usaram arrays de metilação do TCGA para selecionar as regiões do painel.

A consequência é uma **mudança do problema**: em vez de *detectar câncer no sangue*, o experimento passa a *classificar o tecido*. Essa diferença é tratada como limitação e explorada na discussão (seção 6).

---

## 4. Parte III — Experimento com dados de tecido (TCGA)

### 4.1 Dados e amostragem

- **Fonte:** UCSC Xena (Goldman et al., 2020), arquivos `TCGA.<COORTE>.sampleMap/HumanMethylation450.gz`, mais o *probe map* `illuminaMethyl450_hg19_GPL16304_TCGAlegacy` (CpG → gene e cromossomo).
- **Medida:** *beta value* de cada CpG, entre 0 e 1, que corresponde à fração das cópias de DNA metiladas naquele sítio.
- **Tipos de câncer:** cólon (COAD), estômago (STAD), pulmão (LUAD), fígado (LIHC), mama (BRCA) e rim (KIRC). São tipos frequentes, de órgãos diferentes, com tecido normal disponível. Incluímos dois órgãos do trato gastrointestinal (cólon e estômago) como par potencialmente difícil.
- **Rótulos:** derivados do código da amostra `TCGA-XX-YYYY-TT`. Os 12 primeiros caracteres identificam o **paciente**. `TT = 01` indica tumor primário e `TT = 11`, tecido normal adjacente.
- **Amostragem:** incluímos todas as amostras normais, os tumores pareados desses normais (mesmo paciente) e tumores aleatórios até completar 150 por coorte (semente fixa). Isso limita o desbalanceamento entre tipos (a coorte de mama tem 790 tumores) e mantém a matriz dentro da memória do Colab gratuito.

| Tecido | Normal | Tumor | Total |
|---|---|---|---|
| Cólon | 38 | 150 | 188 |
| Estômago | 2 | 150 | 152 |
| Fígado | 50 | 150 | 200 |
| Mama | 98 | 150 | 248 |
| Pulmão | 32 | 150 | 182 |
| Rim | 160 | 160 | 320 |
| **Total** | **380** | **910** | **1.290** |

As 1.290 amostras vêm de **919 pacientes**, dos quais **371 têm tumor e normal**. O rim tem 160 tumores porque todos os seus normais são pareados.

### 4.2 Pré-processamento

1. **Só CpGs autossômicos** (385.809 de 395.985). A coorte de mama é quase só feminina, e os CpGs dos cromossomos sexuais permitiriam ao modelo aprender sexo em vez de tecido.
2. **Remoção de CpGs com algum valor ausente** (49.549), restando 336.260.
3. **Filtro não supervisionado:** os 20 mil CpGs de maior variância. Esse filtro não usa rótulos. Mesmo assim, a auditoria (seção 5.5) mostra que aplicá-lo dentro de cada *fold* de treino produz o mesmo resultado.

Na **análise exploratória por PCA** (PC1 com 34,7% e PC2 com 11,6% da variância), as amostras **normais formam grupos compactos por tecido**, e os **tumores se dispersam**. Isso reflete a heterogeneidade e a desregulação da metilação no câncer. O sinal de tecido é visível antes de qualquer modelo supervisionado.

`[Figura: tcga_pca.png | PCA dos 20 mil CpGs de maior variância, colorida por tecido e tipo (esquerda) e por tumor/normal (direita).]`

### 4.3 Classificador

**Algoritmo: Random Forest** (Breiman, 2001). É um conjunto de árvores de decisão, cada uma treinada numa amostra *bootstrap* dos dados e considerando, a cada divisão, um subconjunto aleatório de features. A predição final é a votação das árvores. Foi escolhido porque:
- lida bem com alta dimensionalidade (p ≫ n);
- captura interações e relações não lineares;
- é pouco sensível à escala das features;
- fornece medidas de importância das features.

Usamos 500 árvores e `class_weight="balanced"`.

**Pipeline** (`scikit-learn`; Pedregosa et al., 2011):
1. **Seleção supervisionada:** `SelectKBest` com teste F (ANOVA), mantendo os *k* CpGs mais associados ao rótulo. Por ser supervisionada, fica **dentro do pipeline**: é reajustada apenas com o treino de cada partição (Ambroise & McLachlan, 2002).
2. **Random Forest.**

### 4.4 Validação

- **Validação cruzada aninhada** (Varma & Simon, 2006):
  - o laço externo (5 *folds*) estima o desempenho, e cada amostra é prevista uma única vez, por um modelo que nunca a viu;
  - o laço interno (3 *folds*, só dentro do treino externo) escolhe *k* ∈ {100, 500, 2000}.
- **Agrupamento por paciente** (`StratifiedGroupKFold`, grupo = paciente), nos dois laços. Tumor e normal do mesmo paciente **nunca** ficam um no treino e outro no teste. Sem isso, o modelo poderia acertar reconhecendo a metilação constitucional do indivíduo, e não o estado tumoral. Verificamos que nenhum paciente aparece em treino e teste de um mesmo *fold*.

**Tarefas:**
1. **Etapa 1 — tumor vs normal:** todas as 1.290 amostras. A métrica de seleção interna é a ROC-AUC.
2. **Etapa 2 — tecido de origem:** os 910 tumores, em 6 classes. A métrica de seleção interna é o F1 macro.
3. **Fluxo hierárquico:** combina as predições fora da amostra das duas etapas, como no artigo. Para os normais erroneamente chamados de tumor, o tecido é previsto por um modelo treinado com todos os tumores, que não inclui nenhum normal. O resultado é uma classificação em 7 classes (normal + 6 tecidos).

**Métricas:** matriz de confusão, precisão, sensibilidade (*recall*), especificidade, F1-score, acurácia e ROC-AUC. No caso multiclasse, as métricas são calculadas por classe (um-contra-todos), e a AUC é a média macro.

**Importância das features:** *permutation importance* num conjunto de teste agrupado por paciente (25% dos pacientes), com o modelo ajustado nos 75% restantes (k = 500). Mede a queda de acurácia ao embaralhar cada CpG. Os CpGs são anotados com o gene mais próximo pelo *probe map*.

---

## 5. Resultados

### 5.1 Etapa 1 — tumor vs normal

Acurácia de **0,984** e ROC-AUC de **0,997**. Em 4 dos 5 *folds* externos, o laço interno escolheu k = 2000.

| Classe | n | Precisão | Sensibilidade | Especificidade | F1 |
|---|---|---|---|---|---|
| Normal | 380 | 0,969 | 0,979 | 0,987 | 0,974 |
| Tumor | 910 | 0,991 | 0,987 | 0,979 | 0,989 |

Desempenho dentro de cada tecido:

| Tecido | Sensibilidade | Especificidade |
|---|---|---|
| Cólon | 1,000 | 0,921 |
| Estômago | 0,993 | 0,000 (apenas 2 normais) |
| Fígado | 0,967 | 0,960 |
| Mama | 0,993 | 0,990 |
| Pulmão | 0,987 | 1,000 |
| Rim | 0,981 | 1,000 |

A sensibilidade é de pelo menos 0,97 em todos os tecidos. A especificidade do estômago não tem valor estatístico, porque o TCGA tem apenas duas amostras normais gástricas.

`[Figura: tcga_etapa1.png | Etapa 1 (tumor vs normal): matriz de confusão e curva ROC das predições da validação cruzada aninhada.]`

### 5.2 Etapa 2 — tecido de origem

Acurácia de **0,980**, F1 macro de **0,980** e ROC-AUC macro (um-contra-todos) de **0,999**.

| Tecido | Precisão | Sensibilidade | Especificidade | F1 |
|---|---|---|---|---|
| Cólon | 0,979 | 0,940 | 0,996 | 0,959 |
| Estômago | 0,929 | 0,960 | 0,986 | 0,944 |
| Fígado | 1,000 | 1,000 | 1,000 | 1,000 |
| Mama | 0,987 | 1,000 | 0,997 | 0,993 |
| Pulmão | 0,993 | 0,987 | 0,999 | 0,990 |
| Rim | 0,994 | 0,994 | 0,999 | 0,994 |

Fígado, mama, pulmão e rim ficam acima de 0,99. Os erros se concentram no par **cólon ↔ estômago**.

`[Figura: tcga_etapa2.png | Etapa 2 (tecido de origem): matriz de confusão das predições da validação cruzada aninhada.]`

### 5.3 Fluxo hierárquico

Acurácia de **0,971** em 7 classes.
- **Normais:** 372 de 380 foram classificados corretamente. Os 8 restantes foram chamados de tumor de cólon (3), estômago (2), fígado (2) e mama (1).
- **Tumores chamados de normal:** 12 (fígado 5, rim 3, pulmão 2, estômago 1, mama 1).
- **Maior confusão:** 9 tumores de cólon foram chamados de estômago, e 3 de estômago foram chamados de cólon.
- **Por classe:** o F1 vai de 0,938 (estômago) a 0,990 (mama).

`[Figura: tcga_hierarquico.png | Fluxo hierárquico (normal + 6 tecidos): matriz de confusão out-of-fold.]`

### 5.4 CpGs mais relevantes

**Etapa 2 (tecido de origem).** Os CpGs de maior importância incluem sítios associados a KANK1 (dois CpGs), ERI3, RBP5, LITAF e TCF3. O beta médio por tecido mostra padrões tecido-específicos:
- **KANK1:** metilado no cólon (0,67–0,74), intermediário no estômago (0,23–0,27) e desmetilado nos demais (≤ 0,16). Funciona como marcador de cólon.
- **RBP5 e LITAF:** desmetilados no fígado (0,31–0,39) e metilados nos demais tecidos. Funcionam como marcadores de fígado.
- **ERI3 e TCF3:** baixos no cólon e altos nos demais. No ERI3, o estômago fica em posição intermediária (0,46).

O **estômago ocupa uma posição intermediária entre o cólon e os demais tecidos** nesses marcadores. Isso é coerente com a confusão cólon ↔ estômago e com a origem embrionária comum dos epitélios gastrointestinais (endoderma).

**Etapa 1 (tumor vs normal).** Destacam-se CpGs próximos de B3GNTL1, FAM78B, HMGA2, FMNL2, MIR646HG e COL8A1.

Os valores de *permutation importance* são baixos (≤ 0,007) porque há **muitos CpGs redundantes**: retirar um único CpG não altera o modelo. Isso indica que o sinal está distribuído por muitos sítios, como confirma a auditoria a seguir.

`[Figura: tcga_heatmap_top_etapa2.png | Beta médio, por tecido, dos 10 CpGs mais importantes da etapa 2 (permutation importance).]`

### 5.5 Auditoria de vazamento

Desempenhos próximos de 100% exigem verificação. Cada teste abaixo usa validação cruzada agrupada por paciente (5 *folds*, k = 500).

| Teste | Etapa 1 | Etapa 2 | Interpretação |
|---|---|---|---|
| Referência (mesmo pipeline) | acc 0,982 · AUC 0,997 | acc 0,978 · AUC 0,999 | — |
| IDs duplicados / pacientes em treino e teste | 0 / 0 | 0 / 0 | Separação por indivíduo correta |
| Filtro de variância ajustado só no treino | acc 0,982 | acc 0,977 | O filtro prévio não inflou o resultado |
| 20 mil CpGs aleatórios | acc 0,988 | acc 0,980 | O sinal está espalhado pelo genoma |
| Rótulos embaralhados | AUC 0,504 | AUC 0,503 (acc 0,18 ≈ 1/6) | O pipeline não fabrica sinal |
| Árvore de decisão com 3 níveis | acc 0,970 | acc 0,647 | Etapa 1 quase trivial; etapa 2 exige muitos CpGs |
| Um único CpG (árvore de 1 nível) | acc 0,905 | acc 0,331 | Um CpG separa ~90% de tumor vs normal |

A auditoria **não encontrou vazamento**. O desempenho elevado reflete uma tarefa genuinamente fácil em tecido:
- a maior parte do DNA de uma amostra tumoral vem das células tumorais;
- cada órgão mantém um programa de metilação próprio, que persiste no tumor. A célula de origem domina a classificação molecular dos tumores (Hoadley et al., 2018).

Um risco que não pode ser testado com os dados disponíveis é o **efeito de lote/placa de processamento** do TCGA associado ao tipo de câncer, porque os identificadores do Xena não trazem a placa. Esse risco é mitigado pela coerência biológica dos resultados (normais agrupados por tecido, marcadores tecido-específicos) e é registrado como limitação.

---

## 6. Parte IV — Discussão

As respostas abaixo seguem a numeração das questões do enunciado.

**Questão 11. Há diferenças de metilação capazes de separar as classes?**
- **Em tecido, sim, e de forma muito robusta.** Tumor e normal são separados com AUC de 0,997, e os seis tecidos de origem com acurácia de 0,98. 20 mil CpGs escolhidos ao acaso produzem o mesmo resultado, e um único CpG já separa 90% das amostras tumor vs normal.
- **No plasma, com o GSE186575, não de forma detectável para câncer vs saudável.** Nenhum promotor passou na correção para testes múltiplos, e o classificador ficou no nível do acaso. Para colorretal vs gástrico, há um sinal multivariado fraco (AUC de 0,79), sem nenhum promotor individualmente significativo.

**Questão 12. Quais regiões contribuem mais?**
- **No tecido**, CpGs com padrão **tecido-específico**: KANK1 (metilado no cólon), RBP5 e LITAF (desmetilados no fígado), ERI3 e TCF3 (desmetilados no cólon). O estômago tem valores intermediários entre o cólon e os demais tecidos, o que explica a principal confusão do modelo.
- **No plasma**, os promotores mais importantes eram majoritariamente pseudogenes ribossomais, o que sugere ruído técnico de mapeamento, e não biologia tumoral.
- Em ambos os casos, a importância individual é baixa por causa da redundância entre features.

**Questão 13. Desempenho e variação entre classes.**
- **No tecido**, fígado, mama, pulmão e rim têm F1 ≥ 0,99 na etapa 2. O cólon tem a menor sensibilidade (0,94), principalmente por confusão com o estômago, que é o par anatômica e embriologicamente mais próximo.
- Na etapa 1, os tumores de fígado são os mais chamados de normal (sensibilidade de 0,967), e a especificidade do estômago não pode ser avaliada (2 normais).
- **No plasma**, a detecção falha (especificidade de 0,32), e a predição de origem é moderada.

**Questão 14. Limitações: tamanho amostral, desbalanceamento e alta dimensionalidade.**
- **Alta dimensionalidade (p ≫ n):** são 20 mil CpGs para 1.290 amostras no TCGA, e 53 mil promotores para 98 amostras no plasma. O risco de sobreajuste foi mitigado pela seleção dentro da validação e pela validação aninhada. No TCGA, a auditoria com rótulos embaralhados confirma que não há sinal espúrio.
- **Desbalanceamento:**
  - o estômago tem apenas 2 normais e a mama tinha 790 tumores;
  - limitamos os tumores a 150 por coorte, usamos `class_weight="balanced"` e reportamos métricas por classe.
- **Tamanho amostral no plasma:** com 31 controles, uma especificidade observada de 100% é compatível com valores reais em torno de 89% (IC de 95%). É impossível verificar os 99% exigidos para rastreamento.
- **Confundimento:**
  - no plasma, o lote e a idade estão associados à classe;
  - no TCGA, não é possível descartar um efeito de lote/placa associado à coorte;
  - os tumores do TCGA são casos já diagnosticados, sem controle de estádio.

**Questão 15. Diferenças em relação à metodologia original.**

| Aspecto | Liu et al. (2020) | Este trabalho — TCGA | Este trabalho — GSE186575 |
|---|---|---|---|
| Material | cfDNA plasmático | Tecido tumoral e normal adjacente | cfDNA plasmático |
| Fração tumoral | Em geral < 1% | Alta (maior parte do DNA) | Em geral < 1% |
| Amostras | 6.689 (multicêntrico, EUA) | 1.290 (919 pacientes) | 98 (um centro, China) |
| Tipos de câncer | > 50 (TOO em ~20 sítios) | 6 | 2 |
| Tecnologia | Bissulfito + captura direcionada | Bissulfito + array 450K | Imunoprecipitação (cfMeDIP) |
| Resolução | CpG por fragmento (~1,1 mi CpGs) | CpG (336 mil, 20 mil usados) | Promotor (~54 mil) |
| Features | *Source models* por região | CpGs selecionados por teste F | Promotores selecionados por teste F |
| Classificador | Regressões logísticas em *ensemble* | Random Forest | Random Forest |
| Limiar | Fixado para especificidade de ~99% | 0,5 / maior probabilidade | 0,5 |
| Validação | Conjunto independente, cego, classificador travado | CV aninhada agrupada por paciente | CV aninhada estratificada |
| Detecção | Especificidade 99,3%; sensibilidade 55% (I–IV) | AUC 0,997 (tumor vs normal) | AUC 0,55 |
| Tecido de origem | Acurácia de 93% | Acurácia de 0,98 (6 classes) | AUC de 0,79 (2 classes) |

A comparação mais importante é entre as colunas: **a assinatura de metilação é forte na fonte (o tumor), mas fica diluída no plasma**. Lá, a fração tumoral é pequena e o sinal é facilmente superado por variação técnica e demográfica. O artigo vence essa distância com três recursos:
- escala de milhares de participantes;
- sequenciamento direcionado profundo;
- leitura do padrão de metilação molécula a molécula, que permite identificar fragmentos tumorais raros em meio ao DNA sanguíneo.

O experimento em tecido também tem valor próprio. Classificadores de metilação em tecido são usados na prática para identificar a origem de **cânceres de sítio primário desconhecido**, que respondem por 3–5% dos cânceres (Moran et al., 2016), e para classificar tumores do sistema nervoso central (Capper et al., 2018). Além disso, o experimento funciona como **controle positivo** do pipeline: mostra que o método funciona quando o sinal existe. Assim, a falha no plasma pode ser atribuída aos dados, e não ao método.

**Questão 16. Por que o resultado não pode ser interpretado como um teste clínico de detecção precoce.**
1. **Material inadequado:** o experimento principal usa **tecido**, cuja obtenção já exige uma biópsia, ou seja, já pressupõe a suspeita de uma lesão. Um teste de rastreamento precisa funcionar no sangue.
2. **Desenho caso-controle:** os pacientes já tinham diagnóstico. Um teste de rastreamento precisa ser avaliado numa população assintomática, onde a prevalência é baixa (cerca de 1,3% ao ano) e o valor preditivo positivo depende fortemente da especificidade.
3. **Especificidade não demonstrada na escala necessária:** com algumas centenas de normais (ou 31 controles, no plasma), não é possível demonstrar uma taxa de falsos positivos abaixo de 1%.
4. **Ausência de validação externa, cega e prospectiva**, e de um classificador travado antes da avaliação. Há ainda possíveis efeitos de lote não testáveis.
5. **Poucos tipos de câncer** e nenhum controle com doenças benignas que alteram a metilação (inflamação, por exemplo).
6. **Sem estratificação por estádio:** a detecção precoce depende do desempenho no estádio I, em que a fração tumoral no plasma é mínima.
7. **Sem demonstração de benefício clínico**, como redução de mortalidade, nem análise de custos e danos.

---

## 7. Conclusão

Este trabalho reproduziu, em pequena escala e com dados públicos, a estrutura de classificação em duas etapas de Liu et al. (2020).

**Tecido (TCGA).** A metilação do DNA separa tumor de tecido normal (AUC de 0,997) e identifica o tecido de origem entre seis tipos de câncer (acurácia de 0,98) com um Random Forest simples. A validação aninhada foi agrupada por paciente, e a seleção supervisionada ficou restrita ao treino. Uma auditoria confirmou que o desempenho não decorre de vazamento: a assinatura de tecido é forte e distribuída por todo o genoma. Os CpGs mais informativos têm padrões tecido-específicos, e a principal confusão (cólon ↔ estômago) é biologicamente coerente.

**Plasma (GSE186575).** Com dados públicos de cfDNA, a detecção de câncer não foi possível: efeitos de lote e idade superaram o sinal tumoral.

**Contraste.** Essa diferença ilustra o desafio central que o artigo enfrenta. Detectar câncer no sangue exige recuperar um sinal raro e diluído, e isso só se torna viável com grande escala amostral, sequenciamento direcionado de alta profundidade, leitura por fragmento e validação independente e cega. Por fim, a demonstração com rótulos aleatórios (AUC de 0,84 com seleção fora da validação, contra 0,27 com seleção dentro) mostra concretamente por que a prevenção de vazamento é indispensável em dados com muito mais features do que amostras.

---

## Registro das fontes, classes e validação

| Item | Experimento principal — TCGA (tecido) | Avaliação experimental — GSE186575 (plasma) |
|---|---|---|
| Dados de metilação | UCSC Xena: `TCGA.{COAD,STAD,LUAD,LIHC,BRCA,KIRC}.sampleMap/HumanMethylation450.gz` | GEO: `GSE186575_DNA-Met_CPM_matrix_promoter.txt.gz` e `GSE186575_series_matrix.txt.gz` (FTP do NCBI) |
| Anotação | *Probe map* `illuminaMethyl450_hg19_GPL16304_TCGAlegacy` (CpG → gene, cromossomo) | Ensembl BioMart (GRCh38): gene, cromossomo, biotipo |
| Classes | Etapa 1: tumor primário (`01`, n = 910) vs normal adjacente (`11`, n = 380). Etapa 2: cólon, estômago, pulmão, fígado, mama, rim (150–160 tumores cada) | Detecção: câncer (CRC 32 + gástrico 35) vs saudável (31). Origem: CRC vs gástrico |
| Validação | CV aninhada com `StratifiedGroupKFold` agrupado por paciente (externa 5, interna 3); seleção de features dentro do pipeline; auditoria de vazamento | CV aninhada estratificada (externa 5 × 3 repetições, interna 5); seleção dentro do pipeline; uma amostra por indivíduo |
| Reprodução | `notebooks/trabalho_final.ipynb`, seção 3 | `notebooks/trabalho_final.ipynb`, seção 4 |

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
9. Hoadley, K. A. et al. (2018). Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell*, 173(2), 291–304.
10. Moran, S. et al. (2016). Epigenetic profiling to classify cancer of unknown primary: a multicentre, retrospective analysis. *Lancet Oncology*, 17(10), 1386–1395.
11. Capper, D. et al. (2018). DNA methylation-based classification of central nervous system tumours. *Nature*, 555, 469–474.
12. Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32.
13. Varma, S. & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, 7, 91.
14. Ambroise, C. & McLachlan, G. J. (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS*, 99(10), 6562–6566.
15. Leek, J. T. et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*, 11, 733–739.
16. Johnson, W. E., Li, C. & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118–127.
17. Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825–2830.
