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
