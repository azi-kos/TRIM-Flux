# Laboratorijski dnevnik — TRIM-Flux

Diplomska naloga, FKKT/FRI, Univerza v Ljubljani  
Mentor: prof. Miha Moškon  
Avtor: Aljaž Koš

---

## 28.4.2026 — Vzpostavitev repozitorija

### Kaj sem naredil
Naredil fork originalnega TRIM repozitorija (uhlerlab/TRIM) na GitHub pod računom azi-kos/TRIM-Flux. Repozitorij sem klonal lokalno in nastavil dva remote-a: `origin` (moj fork) in `upstream` (original). Ustvaril sem 4 branche: `main`, `replicate-original`, `flux-extension`, `rna-flux`.

### Opažanja in rezultati
- Repozitorij vsebuje originalno kodo: `trim.py` (glavni VAE model), `learn_tcr_embedding.py` (CNN za TCR), `analysis/` (evaluacijske skripte)
- `.mm` fajlov ki jih originalna koda bere ni na GEO — avtorji so objavili `data_processing.py` kot primer kako iz surovih podatkov narediti vhodne `.pkl` fajle
- GEO dataset GSE200996 vsebuje: `.h5` fajle (RNA count matrike), `filtered_contig_*.csv.gz` (TCR sekvence), `*meta.data.txt` (metadata po celici)

---

## 1.5.2026 — Preprocessing podatkov

### Kaj sem naredil
Napisal in dokončal preprocessing notebook (`notebooks/01_preprocessing.ipynb`) ki iz surovih GEO podatkov naredi `.pkl` fajle za TRIM. Pipeline naloži RNA `.h5` fajle z `scanpy`, poveže barcodes z metadata, filtrira T celice, naloži TCR kontige, zgradi indekse in shrani vse potrebne fajle.

### Opažanja in rezultati
- **Barcode kolizija**: isti 16-bazni barcode se pojavi pri več pacientih (ker so multipleksirani). Rešitev: povsod uporabljamo kompozitni ključ `(barcode_clean, Patient_ID)` namesto samega barcoda.
- bili so vmesni probelmi z RAM crashi, vendar je bilo to popravljeno
- **TCR dedupliciranje**: `drop_duplicates` je moral priti šele po združitvi vseh fajlov in po dodelitvi `Patient_ID`, sicer so se barcodi različnih pacientov pomešali. Globalni sort po UMI pred `drop_duplicates` zagotovi da obdržimo sekvenco z največ dokazi.
- **Originalni pipeline ima 8 stolpcev** v `data_labels`, ne 4: `Tissue, Treatment Stage, SubCellType, Patient, CDR3(Beta1), tcr_v, tcr_j, treatment`. Vse `analysis/` skripte in oba modela jih unpacked po poziciji. Dodali smo manjkajoče stolpce: `SubCellType` (CD4=0, CD8=1 iz `CellClass`), `tcr_v`/`tcr_j` (indeksi V in J genov iz TCR kontigov), `treatment` (0 — v originalu nikoli definiran, nikjer se ne uporablja).

### Rezultat
- `data_rna.pkl`: sparse matrika `(146776, 36601)`, log1p normalizirana
- `data_labels.pkl`: DataFrame `(146776, 8)` — vsi stolpci v pravilnem vrstnem redu
- `df_all_tcrs.pkl`: `48103` unikatnih CDR3β sekvenc, paddane na enako dolžino

---

## 1.5.2026 — TCR embedding (CNN avtoencoder)

### Kaj sem naredil
Prepisal `learn_tcr_embedding.py` v notebook (`notebooks/02_tcr_embedding.ipynb`) in pognal na Colabu A100. Zamenjal `PCA` z `TruncatedSVD` (sparse matrika).

### Rezultat
- 50100 korakov, loss konvergiral: ~0.01 → ~0.0001
- `data_tcr.pkl`: `(146776, 100)` — TCR embeddingi za vsako celico
- UMAP ne kaže ločevanja po tkivu/fazi — pričakovano, ker CNN trenira samo na rekonstrukciji sekvence, ne na biološkem kontekstu. Ločevanje bo naloga TRIM VAE.

---
## 13.5.2026 — TRIM VAE trening (pacient 0)

### Kaj sem naredil
Prepisal `trim.py` v notebook (`notebooks/03_trim.ipynb`) s prilagoditvami za Colab A100 (argparse → spremenljivke, `cuda:1` → `cuda:0`, `PCA` → `TruncatedSVD`, `os.mkdir` → `os.makedirs`). Pognal trening za heldout pacienta 0 (leave-one-out).

### Težave
- **Pseudo-kloni (pairwise L1 razdalje na 146k × 146k)**: preveč RAM in časa. Rešeno z vzorcem 100k celic — threshold kalibriran na vzorcu, vrednost (`thresh_fitted=1.89`) shranjena v `args.txt`.
- **UMAP transform** na vseh 146k napovedih precount — rešeno z `preds_rna[mask_leave_out]` (samo heldout celice).
- `preds.npz` je bil prvič shranjen prazen (shranjevanje zakomentirano) — ponovljen zagon z naloženim `model.pth`.

### Rezultat
- 20000 korakov, loss konvergiral: `1.092 → 0.108`
- `model.pth` shranjen pri korakih 10000 in 20000
- `preds.npz`: `preds_rna` in `preds_tcr` za heldout pacienta 0
- UMAP vizualizacija v teku

---

## 20.6.2026 — Diagnoza prvih rezultatov in popravek reducerja (veja `replicate-original`)

### Diagnoza
Pregled `umap_rna_patient0.png` je pokazal prazna Tumor panela in ne-prekrivanje napovedi z realnimi celicami. Analiza je razkrila:
- **Pacient 0 nima tumorskih celic** (Blood-Pre=571, Blood-Post=449, Tumor=0) — prazna panela nista bug, ampak napačna izbira heldouta. Le 4 pacienti (P18, P19, P24, P27) imajo vse 4 kvadrante; 22 ima Tumor-Post. Najboljši heldout: P24 (8152 Tumor-Post celic).
- **Trening je zdrav** — loss lepo konvergira (0.285 → 0.108), RNA rekonstrukcija pa je najšibkejši člen (0.329 → 0.218), kar kaže na necentriran SVD prostor.
- **Reducer neujemanje**: notebooka sta uporabljala `TruncatedSVD` (necentriran), original (`trim.py`, `analysis/`) pa `sklearn.PCA` (centriran). Napovedi so zato v drugem prostoru kot pri originalu; pri evalvaciji z `analysis/2.2.eval_gen.py` (ki naredi svoj PCA) bi to dalo nesmiselne rezultate.

### Kaj sem naredil (P1)
- `02_tcr_embedding.ipynb` in `03_trim.ipynb`: zamenjal `TruncatedSVD` → `sklearn.PCA(100)` na gosti `float32` matriki (sparse → dense, `del` + `gc` za RAM). Vrh ~40-50 GB, varno pri High-RAM.
- Shranjevanje **PCA objekta** (`rna_pca.pkl`), da evalvacija dela v istem prostoru (`inverse_transform`).
- Opozorilo: pred ponovnim zagonom je treba **izbrisati star SVD cache** (`data_rna_pca.pkl`, `umap_trained_rna.pkl`) na Drive.

### Kaj sem naredil (P4)
- `03_trim.ipynb`: `heldout_patient = [24]` namesto `[0]`. P24 ima vse 4 kvadrante (Blood-Pre 1548, Blood-Post 2109, Tumor-Pre 5568, Tumor-Post 8152) — smiseln za evalvacijo. Rezultati gredo v `holdout24`.

### Opomba: P3 odpade
Prvotno načrtovan P3 (shrani `tcr_dists` v `preds.npz`) je nesmiseln: `tcr_dists` je `pairwise_distances(preds_tcr)` čez vseh ~146k celic = 146776² × 8B ≈ **172 GB**. Original ga zato v `trim.py` **ne shrani** (zakomentiran v `np.savez`). `2.2.eval_gen.py:184` ga sicer bere, a je to nekonsistentna/mrtva koda v originalu. `our_expansion_prediction` (helpers.py) razdalje računa interno na majhnih generiranih podmnožicah. Notebook 03 torej že sledi originalu — `tcr_dists` se NE shranjuje.

### Kaj sem naredil (P2)
- `01_preprocessing.ipynb`: `df_all_tcrs` dobi 4 count stolpce (`'1'..'4'` = klonska velikost po kondiciji: blood-pre, blood-post, tumor-pre, tumor-post), kot original (`data_processing.py:115-135`). Counte izpeljem iz `data_labels` (CDR3(Beta1) × Tissue × Treatment Stage) in jih **poravnam na obstoječi vrstni red** `df_all_tcrs.index` (NE sortiram — sicer se TCR indeksi razsujejo). NaN za kondicije brez pojavitve (eval dela `.fillna(0)`).
- Preverjeno na dejanskih podatkih: `(48110, 4)`, vsota countov = 146776 = št. celic. Eval-stil dostopi (`iloc[:,1]>iloc[:,0]`, `.sum(axis=1)`) delajo.
- **Brez tega so njihove metrike ekspanzije (clont_count, expansion ROC) crashale.** Zdaj originalne `analysis/` metode tečejo skoraj nespremenjene.

### Kaj sem naredil (pseudo-kloni — crash fix)
- `03_trim.ipynb` cell-25: `tcr_dists` se je računal na **vseh ~146k celicah** (`pairwise_distances(preds_tcr, preds_tcr)`) = 146776² × 8B ≈ **172 GB → crash** (isto kot 13.5.). Original to računa na 2 TB stroju; Colab High-RAM ima le ~83 GB system RAM (A100 80GB je VRAM, ne system RAM), zato večji Colab to ne reši.
- **Popravek:** `tcr_dists` + `get_pseudoclones` zdaj delata **samo na heldout celicah** (P24 = 17377 celic → matrika 2.42 GB). Ker `get_pseudoclones` maskira `dists[i]` na isti `patient`+`bloodtumor`, je rezultat za heldout celice **identičen** polni matriki, brez odpadka. `pseudo_tcrs` se razširi nazaj na polno dolžino (-10 placeholder), da `preds.npz[mask_leave_out]` deluje.
- Preverjeno: P24 ima 1548 blood-pre celic, 497 unikatnih klonov za kalibracijo `thresh` (dovolj), vse 4 kondicije prisotne.

### Naslednji korak
- P5: pognati originalne evalvacijske metrike iz `analysis/HNSCC/` (clont_count_by_cluster Pearson r, diversity r) na novih `preds.npz` (heldout P24). Treba: prilagoditi vhode (`.pkl`, `ID_NULL_TCR=-1`, en heldout namesto 27, PCA cache).
- **Pomembno:** `df_all_tcrs` zdaj ima stolpce — pred ponovnim zagonom notebooka 02/03 ni potrebnih sprememb (oba bereta le `.index`), a `df_all_tcrs.pkl` na Drive je treba **regenerirati z notebookom 01**.

---

## 20.6.2026 — Notebook 04: evalvacija klonske ekspanzije (Figure 4)

### Kaj sem naredil
Ustvaril `notebooks/04_expansion_eval.ipynb` — evalvacija napovedi klonske ekspanzije za **enega heldout pacienta naenkrat** (P24). Logika **prekopirana iz originala** (`analysis/HNSCC/helpers.py`), prilagojena na študentove vhode:
- **`baseline_expansion_prediction_single`** (iz `helpers.py:284`): navaden klasifikator (kNN/NN × rna/tcr/both) na realnih blood celicah → ROC AUC. Ne rabi modela.
- **`our_expansion_prediction_single`** (iz `helpers.py:395`): iz blood-pre `z` heldout pacienta generira celice v 4 kondicijah, šteje pseudo-klonalnost (`thresh_fitted/2`), primerja z resnično (`df_all_tcrs` counti).
- **Končna ROC AUC** (iz `2.2.eval_gen.py:569-575`): `x = blood-post > blood-pre` (resnica), `y = pseudo BA/BB` (napoved).

### Premoščena neujemanja (original → študent)
`.npz`→`.pkl`, `data_labels` DataFrame→`.values`, `ID_NULL_TCR` 0→-1, trde poti→Drive, globalni scope, `Generator()` brez args (študentova verzija), **izpuščen mrtvi `tcr_dists`** (helpers.py:446).

### Validacija (lokalno, s simuliranimi P2 counti)
- JSON + sintaksa vseh celic: OK
- Baseline P24: train 56795 celic, test 1503 (po OOD filtru), **60 ekspandiranih → ROC definiran**
- Our P24: 1548 blood-pre, `real_clon` (1548,4), **69 ekspandiranih → ROC definiran**
- Poravnava `recon_*_z` (03 `mask_leave_out`) ↔ `x_label_ho` (04): identična ✅

### Pogoji za zagon (Colab)
- `holdout24/preds.npz` z `recon_*_z` (cell-27 mora teči po pseudo-klonih)
- `df_all_tcrs.pkl` na Drive = **P2-verzija** (regeneriraj z notebookom 01 — lokalna kopija je še stara, brez countov)
- Per-patient: ROC na enem pacientu = šibka statistika; agregat čez 3-5 pacientov kasneje (razširi `output_folders`)

---

## 27.6.2026 — Flux razširitev: raziskava izvedljivosti + scFEA demo (veja `flux-extension`)

### Merge
Vse replikacijsko delo (P1, P2, P4, pseudo-fix, notebook 04) mergano v `main` (fast-forward). Vse veje (`replicate-original`, `flux-extension`, `rna-flux`) poravnane na `fe529ff`. Flux delo se nadaljuje na `flux-extension`.

### Raziskava izvedljivosti flux pipeline
Cilj Var 2 (RNA+TCR+Flux). Vhod: 146776 T-celic. Raziskal orodja (CLAUDE.md je predlagal iMAT/ftINIT — suboptimalno):
- **Compass** (FBA): ~30 min/vzorec → per-cell na 147k = leta. NE sprejme zunanjega GEM. Izključen.
- **ftINIT** (RAVEN/MATLAB): gradi context-specific MODEL, ne per-cell flux. Ni naš pipeline.
- **scFEA** (GNN, GPU): IZBRAN. Per-cell, GPU-pospešen, ~168 modulov za človeka. Dropout-robusten (korelacija >0.85 pri simuliranem dropoutu, testirano v članku). Ne klasificira celičnega tipa → ne meša celic. Brez poolanja.

### Razčiščene skrbi (vse iz virov)
- Determinizem: znotraj zagona da; fiksiraj seed. Stabilnejši od Compass (LP ima več optimumov).
- Mešanje celic: nemogoče — flux iz lastne ekspresije, ne klasifikacije.
- ftINIT→Compass ideja odpadla: Compass ne sprejme zunanjega modela.

### Kaj sem naredil
Ustvaril `notebooks_flux/01_scfea_demo.ipynb` — proof-of-concept: namesti scFEA, vzame 2000 celic (P24), preveri orientacijo (geni×celice) + gene imena (simboli vs Ensembl), požene scFEA, **izmeri čas** → ekstrapolacija na cel dataset, pregleda flux izhod (celice × ~168 modulov).

### Odprto / pred zagonom
- Preveri **gene imena** v `data_rna` (scFEA rabi simbole, ne Ensembl ID); če Ensembl → pretvorba (mygene/pybiomart).
- Preveri točna imena `module_gene` + `cmMat` v scFEA `data/` (auto-detekcija v cell-11).
- Pridržek za diplomo: scFEA fluksi so RELATIVNI, model-odvisni; benchmark scFEA/Compass/METAFlux ne obstaja.

---