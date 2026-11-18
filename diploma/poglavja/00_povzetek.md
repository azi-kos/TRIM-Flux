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

Imunoterapija z zaviralci imunskih kontrolnih točk je učinkovita le pri delu bolnikov z rakom. Odziv nanjo spremlja klonska ekspanzija limfocitov T, zato bi napoved, kateri kloni se bodo po zdravljenju namnožili, pomagala razumeti in napovedati odziv. Model TRIM, pogojni variacijski avtokodirnik, ekspanzijo napove iz transkriptoma in receptorja celic T (TCR) posameznih celic pred zdravljenjem. Namen dela je bil ugotoviti, ali metabolno stanje celic, ocenjeno iz podatkov sekvenciranja RNA posameznih celic, napovedi doda novo informacijo. Model TRIM smo ponovili na javnih podatkih 28 bolnikov s ploščatoceličnim karcinomom glave in vratu, zdravljenih z neoadjuvantno imunoterapijo (69 930 celic T po kontroli kakovosti), z lastno predobdelavo. Razširili smo ga z metabolnimi fluksi skozi 168 metabolnih modulov, ki smo jih ocenili z metodo scFEA. Tri različice modela (RNA + TCR, fluksi + TCR ter RNA + fluksi + TCR) smo ovrednotili s shemo »izpusti enega bolnika« pri dveh bolnikih, vsako s tremi ponovitvami učenja z različnimi semeni. Različica RNA + TCR je ekspanzijo klonov v krvi napovedala s povprečno ploščino pod krivuljo ROC (AUC) 0,87 in 0,84 (posamezni zagoni od 0,77 do 0,92) in v štirih od petih zagonov presegla izhodiščne klasifikatorje. Velikost klona pred zdravljenjem, ki jo model lahko razbere iz TCR, pa je ekspanzijo napovedala z AUC 0,96 in 0,95, torej bolje od vseh posameznih zagonov modela in enako dobro kot združene napovedi zagonov z različnimi semeni. Model je velikosti klona nekaj dodal le pri klonih z eno celico pred zdravljenjem pri enem bolniku (10 ekspandiranih klonov), in to samo v različicah s transkriptomom. Fluksi namesto transkriptoma so dali slabšo in manj stabilno napoved, fluksi poleg transkriptoma pa je niso zanesljivo izboljšali. Razlike med zagoni z različnimi semeni so bile večje od razlik med različicama s transkriptomom. Sklepamo, da pri uporabljeni oznaki ekspanzije AUC na celicah v veliki meri odraža velikost klona, fluksi scFEA pa v naši zasnovi napovedi niso izboljšali. Za oceno prispevka metabolizma bo treba model ovrednotiti pri več bolnikih, z izhodiščem iz velikosti klona in z merili na klonih z eno celico.

**Ključne besede:** limfociti T, klonska ekspanzija, sekvenciranje posameznih celic, metabolni fluks, variacijski avtokodirnik

::: {custom-style="Naslov povzetka"}
Predicting T-cell clonal expansion by integrating transcriptome, T-cell receptors and metabolic fluxes
:::

**Abstract**

Immune checkpoint inhibitor therapy is effective in only a subset of cancer patients. Response is accompanied by clonal expansion of T cells, so predicting which clones will expand after treatment could help to understand and predict the response. The TRIM model, a conditional variational autoencoder, predicts expansion from the transcriptome and the T-cell receptor (TCR) of single cells before treatment. The aim of this work was to determine whether the metabolic state of cells, estimated from single-cell RNA sequencing data, adds new information to the prediction. We reproduced the TRIM model on public data from 28 patients with head and neck squamous cell carcinoma treated with neoadjuvant immunotherapy (69,930 T cells after quality control), using our own preprocessing. We extended it with metabolic fluxes through 168 metabolic modules, estimated with the scFEA method. Three model variants (RNA + TCR, fluxes + TCR, and RNA + fluxes + TCR) were evaluated with a leave-one-patient-out scheme in two patients, each trained three times with different random seeds. The RNA + TCR variant predicted clonal expansion in blood with a mean area under the ROC curve (AUC) of 0.87 and 0.84 (individual runs 0.77 to 0.92) and outperformed the baseline classifiers in four of five runs. However, clone size before treatment, which the model can infer from the TCR, predicted expansion with an AUC of 0.96 and 0.95, which is better than every individual model run and as good as the combined predictions of runs with different seeds. The model added some information beyond clone size only for clones with a single cell before treatment in one patient (10 expanded clones), and only in the variants that included the transcriptome. Fluxes instead of the transcriptome gave a weaker and less stable prediction, and fluxes in addition to the transcriptome did not reliably improve it. Differences between runs with different seeds were larger than differences between the two variants that included the transcriptome. We conclude that, with the expansion label used, cell-level AUC largely reflects clone size, and that scFEA fluxes did not improve the prediction in our setup. Assessing the contribution of metabolism will require evaluation in more patients, with a clone-size baseline and with metrics on single-cell clones.


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
| GEM | genomski metabolni model (angl. *genome-scale metabolic model*) |
| GEO | zbirka podatkov o izražanju genov (angl. *Gene Expression Omnibus*) |
| GPR | gen–protein–reakcija (angl. *gene–protein–reaction*) |
| HPV | humani papilomavirus (angl. *human papillomavirus*) |
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

| Simbol | Pomen |
|---|---|
| $c_i$ | klon celice $i$ |
| $d$ | razsežnost latentnega prostora |
| $\delta$ | najmanjša želena razdalja med vektorji TCR celic različnih klonov (enačba 5) |
| $\varepsilon$ | majhna konstanta pri normalizaciji štetij (enačba 3) |
| $\eta^2$ | delež variance, ki ga pojasni pripadnost skupini |
| $f_{im}$ | fluks skozi metabolni modul $m$ v celici $i$ (enačba 4) |
| $\hat{t}_i$ | rekonstruiran vektor TCR celice $i$ |
| $\lambda$ | utež člena funkcije izgube |
| $\mathcal{L}$ | funkcija izgube |
| $\mu$, $\sigma$ | povprečje in standardni odklon latentne porazdelitve celice |
| $\bar{p}_i$ | povprečna velikost psevdoklona celice $i$ (enačba 6) |
| $\rho$ | Spearmanov koeficient korelacije |
| $S$ | stehiometrična matrika (enačba 2) |
| $\tau$ | prag razdalje za združevanje vektorjev TCR v psevdoklone |
| $v$ | vektor metabolnih fluksov (enačba 2) |
| $x_{ig}$ | število prepisov gena $g$ v celici $i$ (enačba 3) |
| $y_i$ | ocena ekspanzije celice $i$ (enačba 6) |
| $z$ | latentni vektor |
