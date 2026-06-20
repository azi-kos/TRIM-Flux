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

### Naslednji korak
- P5: kvantitativna metrika uspeha (kNN-overlap predicted vs real tumor-post) za heldout P24.

---