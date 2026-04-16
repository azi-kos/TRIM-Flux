# TRIM-Flux — Projekt za diplomsko nalogo

## Kontekst projekta

To je diplomski projekt na FKKT, Univerza v Ljubljani.
Mentor: Prof. Miha Moškon (FRI).

Cilj je replicirati TRIM model (He et al., Nature Communications 2026)
in ga razširiti z metabolnimi flux podatki kot tretjo modaliteto.

## Kaj je TRIM

TRIM (TCR-RNA Integrating Model) je multimodalni conditional VAE ki:
- Integrira scRNA-seq (genska ekspresija) in scTCR-seq (T celični receptor)
- Napoveduje stanje T celic v tumorju PO imunoterapiji
- Na podlagi krvnih T celic pacienta PRED zdravljenjem
- Apliciran na raku glave in vratu (HNSCC), kolorektalnem raku (CRC), pan-cancer

Originalni repozitorij: https://github.com/uhlerlab/TRIM
Originalni članek: https://www.nature.com/articles/s41467-026-70505-0

## Arhitektura TRIM

### Predprocesiranje
- RNA: surovi scRNA-seq → PCA → 100 dim
- TCR: CDR3β sekvenca → Atchleyevi faktorji (5 dim/AA) → CNN avtoencoder
  (3× Conv1d, stride=2, BatchNorm, LeakyReLU → FC → 100 dim embedding)
  CNN se trenira 50.100 korakov, nato dekoder zavrže, encoder zamrznjen

### Kondicionalni vektorji (128 dim vsak)
- Tri vrste: tkivo (kri=0/tumor=1), čas (pred=0/po=1), pacient ID
- Nastanejo prek MLP: RNA PCA (100) → 2048 → 1024 → 512 → 128
- groupby_mean agregacija: povprečje celic z istim labelom
- Skupaj: 3 × 128 = 384 dim
- KLJUČNO: posodabljajo se z backpropom pri vsaki iteraciji skupaj z modelom
- Prisotni pri ENCODERJU in DECODERJU

### TRIM VAE
Encoderja (MLP):
- Vhod: data(100) + cond_vektorji(384) = 484 dim
- Arhitektura: 484 → 2048 → 1024 → 512 → 2048
- Izhod: mu(1024) + logvar(1024) — za RNA in TCR ločeno

Fusion:
- mu = mean(mu_rna, mu_tcr) → 1024 dim
- logvar = mean(logvar_rna, logvar_tcr) → 1024 dim
- z = reparameterize(mu, logvar) → 1024 dim
- Enostavno aritmetično povprečje, brez attention

Dekoderja (MLP):
- Vhod: z(1024) + cond_vektorji(384) = 1408 dim
- Arhitektura: 1408 → 512 → 1024 → 2048 → 100
- RNA decoder → transkripcijsko stanje
- TCR decoder → klonalna ekspanzija

### Loss funkcija (4 komponente)
- λ=15  KL divergenca (regularizacija latentnega prostora)
- λ=1   MSE rekonstrukcija RNA
- λ=1   Kontrastivna loss za TCR klone
         (isti klon: L1 razdalja ↓, različni kloni: margin loss ↑)
- λ=var Regularizacija embedding norm

Trening: 20.000 korakov, batch=5120 celic, Adam optimizer
Evalvacija: leave-one-patient-out (27 pacientov pri HNSCC)

### Counterfactual inference
1. Vzameš krvne celice pacienta (kri, pred zdravljenjem)
2. Enkodiraš z kondionalnimi vektorji [kri, pred]
3. Dobiš z — "identiteta" celice brez konteksta
4. Dekodiraš z ZAMENJАNIMI kondionalnimi vektorji [tumor, post]
5. Rezultat: napoved stanja tumorskih T celic po zdravljenju

## Podatki

Trije neodvisni dataseti:
- HNSCC: GEO GSE200996 (Luoma et al. 2022, Cell)
- CRC:   GEO GSE236581 (Chen et al. 2024, Cancer Cell)
- Pan-cancer: GEO GSE156728 (Zheng et al. 2021, Science)

Format vhodnih podatkov (pickle fajli):
- data_rna.pkl: numpy array (n_cells, n_genes)
- data_labels.pkl: DataFrame z Tissue, Treatment Stage, Patient, CDR3(Beta1)
- data_tcr.pkl: numpy array (n_cells, 100) — TCR embeddingi
- df_all_tcrs.pkl: DataFrame vseh unikatnih TCR sekvenc

## Moje razširitve (cilj diplome)

### Varianta 1: Replikacija originala
Repliciraj TRIM na HNSCC datasetu, primerjaj z objavljenimi rezultati.

### Varianta 2: TRIM + Flux
Dodaj flux encoder/decoder kot tretjo modaliteto:
- Flux matrika iz iMAT/ftINIT na Human1 GEM
- Flux encoder: flux(500-2000) + cond(384) → MLP → mu_flux, logvar_flux
- Fusion: mean(mu_rna, mu_tcr, mu_flux)
- Flux decoder: z + cond → flux rekonstrukcija
- Loss: dodaj MSE_flux komponento
- Hipoteza: metabolni fluxe dodajo informacijo o exhaustion stanju
  ki je v RNA le posredno prisotna

### Varianta 3: RNA + Flux (brez TCR)
Preveri ali flux nadomesti TCR signal za napoved transkripcijskega stanja.

## Metabolni flux pipeline

Model: Human1 (SysBioChalmers/Human-GEM)
Metoda: iMAT ali ftINIT (RAVEN Toolbox / COBRApy)
Vhod: scRNA-seq count matrika
Izhod: flux matrika (n_cells, n_reactions)

Relevantni metabolni procesi za T celično exhaustion:
- Oksidativna fosforilacija (OXPHOS) ↓ pri exhausted
- Glikoliza — spremenjena regulacija
- Fatty acid oxidation (FAO)
- Arginin/glutamin metabolizem

## Okolje in infrastruktura

- Trening: Google Colab A100 (40GB VRAM)
- Lokalno: VS Code za pisanje kode, Git za verzioniranje
- Python 3.10, PyTorch, COBRApy, scanpy, anndata

## Git workflow

Branche:
- main: stabilna koda
- replicate-original: replikacija TRIM
- flux-extension: dodajanje flux modalitete
- rna-flux: varianta brez TCR