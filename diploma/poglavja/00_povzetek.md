<!--
PRVI ODSEK po Vzorcu 2026: naslovnico, notranjo naslovnico in Izjavo o avtorstvu naredi zgradi.sh samodejno
iz metadata.yaml (predloga/razdelki.py). Zahvala (neobvezna) bi prišla tukaj, pred povzetek.
Vsak naslov s slogom »Naslov povzetka« ali »Naslov brez stevilke« začne novo liho stran in NI v kazalu.
Povzetek: ½ do 1 stran (področje, namen, metodološki pristop, ključni rezultati, kratek zaključek), brez slik in citatov,
razumljiv samostojno, 3–5 ključnih besed. Abstract: največ 1 stran, s prevedenim naslovom, 3–5 keywords.
Napišemo ga ZADNJEGA.
-->

::: {custom-style="Naslov povzetka"}
Napovedovanje klonske ekspanzije limfocitov T z združevanjem transkriptoma, receptorjev T in metabolnih fluksov
:::

**Povzetek**

Imunoterapija z zaviralci imunskih kontrolnih točk je učinkovita le pri delu bolnikov z rakom. Odziv nanjo spremlja klonska ekspanzija limfocitov T, zato bi napoved, kateri kloni se bodo po zdravljenju namnožili, pomagala razumeti in napovedati odziv. Model TRIM, pogojni variacijski avtokodirnik, ekspanzijo napove iz transkriptoma in receptorja celic T (TCR) posameznih celic pred zdravljenjem. Namen dela je bil ugotoviti, ali metabolno stanje celic, ocenjeno iz podatkov sekvenciranja RNA posameznih celic, napovedi doda novo informacijo. Model TRIM smo ponovili na javnih podatkih 28 bolnikov s ploščatoceličnim karcinomom glave in vratu, zdravljenih z neoadjuvantno imunoterapijo (69 930 celic T po kontroli kakovosti), z lastno predobdelavo. Razširili smo ga z metabolnimi fluksi skozi 168 metabolnih modulov, ki smo jih ocenili z metodo scFEA. Tri različice modela (RNA + TCR, fluksi + TCR ter RNA + fluksi + TCR) smo ovrednotili s shemo »izpusti enega bolnika« pri dveh bolnikih, vsako s tremi semeni. Različica RNA + TCR je ekspanzijo klonov v krvi napovedala s ploščino pod krivuljo ROC (AUC) do 0,92 in presegla izhodiščne klasifikatorje, kar je primerljivo z izvirnim delom. Enako dobro (AUC 0,96 in 0,95) pa je ekspanzijo napovedala že velikost klona pred zdravljenjem, ki jo model razbere iz TCR. Model je velikosti klona nekaj dodal le pri klonih z eno celico pred zdravljenjem, in to samo v različicah s transkriptomom. Fluksi namesto transkriptoma so dali slabšo in manj stabilno napoved, fluksi poleg transkriptoma pa je niso zanesljivo izboljšali. Razlike med zagoni z različnimi semeni so bile večje od razlik med različicami. Sklepamo, da pri uporabljeni oznaki ekspanzije AUC na celicah v veliki meri odraža velikost klona, fluksi scFEA pa ne prinesejo informacije, ki je transkriptom ne bi že vseboval. Za oceno prispevka metabolizma bo treba model ovrednotiti pri več bolnikih, z izhodiščem iz velikosti klona in z merili na klonih z eno celico.

**Ključne besede:** limfociti T, klonska ekspanzija, sekvenciranje posameznih celic, metabolni fluks, variacijski avtokodirnik

::: {custom-style="Naslov povzetka"}
Predicting T-cell clonal expansion by integrating transcriptome, T-cell receptors and metabolic fluxes
:::

**Abstract**

Immune checkpoint inhibitor therapy is effective in only a subset of cancer patients. Response is accompanied by clonal expansion of T cells, so predicting which clones will expand after treatment could help to understand and predict the response. The TRIM model, a conditional variational autoencoder, predicts expansion from the transcriptome and the T-cell receptor (TCR) of single cells before treatment. The aim of this work was to determine whether the metabolic state of cells, estimated from single-cell RNA sequencing data, adds new information to the prediction. We reproduced the TRIM model on public data from 28 patients with head and neck squamous cell carcinoma treated with neoadjuvant immunotherapy (69,930 T cells after quality control), using our own preprocessing. We extended it with metabolic fluxes through 168 metabolic modules, estimated with the scFEA method. Three model variants (RNA + TCR, fluxes + TCR, and RNA + fluxes + TCR) were evaluated with a leave-one-patient-out scheme in two patients, each with three random seeds. The RNA + TCR variant predicted clonal expansion in blood with an area under the ROC curve (AUC) of up to 0.92 and outperformed the baseline classifiers, comparable to the original work. However, clone size before treatment, which the model can infer from the TCR, predicted expansion equally well (AUC 0.96 and 0.95). The model added some information beyond clone size only for clones with a single cell before treatment, and only in the variants that included the transcriptome. Fluxes instead of the transcriptome gave a weaker and less stable prediction, and fluxes in addition to the transcriptome did not reliably improve it. Differences between runs with different seeds were larger than differences between variants. We conclude that, with the expansion label used, cell-level AUC largely reflects clone size, and scFEA fluxes do not provide information that the transcriptome does not already contain. Assessing the contribution of metabolism will require evaluation in more patients, with a clone-size baseline and with metrics on single-cell clones.


**Keywords:** T cells, clonal expansion, single-cell sequencing, metabolic flux, variational autoencoder

::: {custom-style="Naslov brez stevilke"}
Kazalo
:::

```{=openxml}
<w:p><w:r><w:fldChar w:fldCharType="begin" w:dirty="true"/></w:r><w:r><w:instrText xml:space="preserve"> TOC \o "1-3" \h \z \u </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>Kazalo se izpolni, ko Word posodobi polja (desni klik → Posodobi polje).</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>
```

::: {custom-style="Naslov brez stevilke"}
Seznam uporabljenih kratic in simbolov
:::

<!-- abecedno, razlaga slovensko (angleško le, če je sicer nejasno); vse kratice iz besedila.
     NE: SI enote, kemijske formule, uveljavljene kratice (DNA, RNA, PCR ...), kratice, uporabljene samo enkrat. -->

| Kratica | Pomen |
|---|---|
| AUC | ploščina pod krivuljo ROC (angl. *area under the curve*) |
| CDR3 | tretja regija, ki določa komplementarnost (angl. *complementarity-determining region 3*) |
| CTLA-4 | antigen 4, povezan s citotoksičnimi limfociti T (angl. *cytotoxic T-lymphocyte-associated protein 4*) |
| FBA | analiza ravnotežja fluksov (angl. *flux balance analysis*) |
| GEM | genomski metabolni model (angl. *genome-scale metabolic model*) |
| GEO | zbirka podatkov o izražanju genov (angl. *Gene Expression Omnibus*) |
| GPR | gen–protein–reakcija (angl. *gene–protein–reaction*) |
| HNSCC | ploščatocelični karcinom glave in vratu (angl. *head and neck squamous cell carcinoma*) |
| kNN | metoda $k$ najbližjih sosedov (angl. *k-nearest neighbours*) |
| NN | nevronska mreža, tu večplastni perceptron (angl. *neural network*) |
| PCA | analiza glavnih komponent (angl. *principal component analysis*) |
| PD-1 | receptor programirane celične smrti 1 (angl. *programmed cell death protein 1*) |
| ROC | krivulja delovanja sprejemnika (angl. *receiver operating characteristic*) |
| scFEA | analiza fluksov v posameznih celicah (angl. *single-cell flux estimation analysis*) |
| scRNA-seq | sekvenciranje RNA posameznih celic (angl. *single-cell RNA sequencing*) |
| scTCR-seq | sekvenciranje receptorjev celic T posameznih celic (angl. *single-cell T-cell receptor sequencing*) |
| TCR | receptor celic T (angl. *T-cell receptor*) |
| TRB | veriga β receptorja celic T (angl. *T-cell receptor beta chain*) |
| TRIM | model za združevanje TCR in RNA (angl. *TCR-RNA Integrating Model*) |
| UMAP | enakomerna aproksimacija in projekcija mnogoterosti (angl. *uniform manifold approximation and projection*) |
| UMI | edinstveni molekulski identifikator (angl. *unique molecular identifier*) |
| VAE | variacijski avtokodirnik (angl. *variational autoencoder*) |
