# Laboratorijski dnevnik – TRIM-Flux

Diplomsko delo, UL FKKT  
Mentor: prof. dr. Miha Moškon  
Avtor: Aljaž Kos

Dnevnik sem 12. 11. 2026 skrajšal. Datumi in vsebina vnosov so ostali enaki, popravki so v zadnjem vnosu. Prejšnja daljša različica je v zgodovini repozitorija (commit 678cb41). Pri programiranju in pisanju dnevnika sem uporabljal orodje UI (Claude).

---

## 28. 4. 2026 – Repozitorij

Naredil sem fork uhlerlab/TRIM (azi-kos/TRIM-Flux), remota `origin` in `upstream` ter veje `main`, `replicate-original`, `flux-extension` in `rna-flux`.

- V repozitoriju so `trim.py` (VAE), `learn_tcr_embedding.py` (CNN za TCR) in `analysis/`.
- Datotek `.mm`, ki jih bere originalna koda, na GEO ni. Avtorji so dali `data_processing.py` kot primer, kako narediti vhodne `.pkl`.
- GSE200996 ima `.h5` (RNA), `filtered_contig_*.csv.gz` (TCR) in `*meta.data.txt`.

## 1. 5. 2026 – Predobdelava

Notebook `01_preprocessing`: nalaganje `.h5` s scanpy, povezava z metapodatki, celice T, kontigi TCR, indeksi, shranjevanje `.pkl`.

- Ista barkoda se pojavi pri več bolnikih (multipleksiranje), zato uporabljam ključ (`barcode_clean`, `Patient_ID`).
- Vmes je zmanjkovalo RAM-a, to sem rešil.
- `drop_duplicates` za TCR mora biti po združitvi datotek in dodelitvi bolnika. Pred tem razvrstim po UMI, da ostane kontig z največ dokazi.
- `data_labels` ima v originalu 8 stolpcev, ne 4. Dodal sem `SubCellType` (CD4 = 0, CD8 = 1), `tcr_v`, `tcr_j` in `treatment` (0, v originalu se ne uporablja).

Rezultat: `data_rna.pkl` (146776 × 36601, log1p), `data_labels.pkl` (146776 × 8), `df_all_tcrs.pkl` (48103 CDR3β).

## 1. 5. 2026 – Vložitev TCR

`learn_tcr_embedding.py` sem prepisal v `02_tcr_embedding` in ga pognal na Colabu (A100). PCA sem zamenjal s `TruncatedSVD` zaradi redke matrike.

- 50100 korakov, izguba z ~0.01 na ~0.0001.
- `data_tcr.pkl`: 146776 × 100.
- UMAP ne loči po tkivu in fazi. To je pričakovano, ker se CNN uči samo rekonstrukcije zaporedja.

## 13. 5. 2026 – Trening TRIM (bolnik 0)

`trim.py` sem prepisal v `03_trim` (argparse v spremenljivke, `cuda:0`, `TruncatedSVD`, `os.makedirs`) in učil z izpuščenim bolnikom 0.

- Matrika razdalj L1 za psevdoklone na 146k celicah ne gre v RAM. Prag sem umeril na vzorcu 100k celic (`thresh_fitted` = 1.89).
- UMAP transform delam le na celicah izpuščenega bolnika.
- `preds.npz` se prvič ni shranil (zakomentirano). Ponovil sem z naloženim `model.pth`.

Rezultat: 20000 korakov, izguba 1.092 do 0.108, `model.pth` pri 10000 in 20000, `preds.npz` za bolnika 0.

## 20. 6. 2026 – Diagnoza prvih rezultatov, PCA namesto SVD

- Bolnik 0 nima tumorskih celic (kri pred 571, kri po 449, tumor 0). Vse štiri skupine imajo le P18, P19, P24 in P27. Za izpuščenega sem izbral P24 (kri pred 1548, kri po 2109, tumor pred 5568, tumor po 8152).
- Izguba konvergira (0.285 do 0.108). Najslabša je rekonstrukcija RNA (0.329 do 0.218).
- Notebooka sta uporabljala `TruncatedSVD` (necentriran), original pa `sklearn.PCA`. V 02 in 03 sem dal `PCA(100)` na gosti float32 matriki (do 40–50 GB RAM) in shranil `rna_pca.pkl`. Star SVD cache na Drive je treba zbrisati.
- `tcr_dists` se ne shranjuje (146776² × 8 B ≈ 172 GB). Tudi original ga ne shrani.
- `01`: `df_all_tcrs` dobi 4 stolpce z velikostjo klona po skupinah, kot v `data_processing.py:115-135`, poravnano na obstoječi indeks. Preverjeno: (48110, 4), vsota 146776. Brez tega so originalne metrike ekspanzije padle.
- `03` cell-25: `tcr_dists` in `get_pseudoclones` računam samo na celicah izpuščenega bolnika (P24: 17377 celic, 2.42 GB). `get_pseudoclones` primerja le celice istega bolnika in tkiva, zato so psevdokloni enaki kot s polno matriko. Colab ima ~83 GB RAM.
- P24 ima 1548 celic kri pred in 497 klonov za umerjanje praga.

Naprej: originalne metrike iz `analysis/HNSCC` na P24, `df_all_tcrs.pkl` znova narediti z 01.

## 20. 6. 2026 – Notebook 04: vrednotenje ekspanzije

`04_expansion_eval` za enega izpuščenega bolnika (P24), logika iz `analysis/HNSCC/helpers.py`:

- `baseline_expansion_prediction_single` (helpers.py:284): kNN in NN na RNA, TCR in obojem.
- `our_expansion_prediction_single` (helpers.py:395): iz z celic kri pred generira celice v štirih skupinah, prešteje psevdoklone (`thresh_fitted/2`) in primerja z resničnimi kloni.
- AUC kot v `2.2.eval_gen.py:569-575`: x = kri po > kri pred, y = razmerje psevdoklonov po/pred.

Prilagoditve: `.npz` v `.pkl`, `ID_NULL_TCR` −1, poti na Drive, brez `tcr_dists`.

Lokalni preizkus s simuliranimi števili: osnovni modeli 1503 testnih celic in 60 ekspandiranih, naš model 1548 celic in 69 ekspandiranih. Poravnava `recon_*_z` z oznakami se ujema. En bolnik je šibka statistika, pozneje več bolnikov.

## 27. 6. 2026 – Fluksi: izbira orodja in preizkus scFEA

Replikacijo sem združil v `main`, vse veje so na `fe529ff`. Fluksi gredo naprej na `flux-extension`.

Orodja za Var 2 (RNA + TCR + fluks) na 146776 celicah. Prvi predlog (CLAUDE.md) je bil iMAT/ftINIT.

- Compass: ~30 min na vzorec, za posamezne celice prepočasen, ne sprejme zunanjega GEM.
- ftINIT (RAVEN, MATLAB): naredi model za kontekst, ne fluksov po celicah.
- scFEA: izbran. Po celicah, GPU, 168 modulov. V članku robusten na dropout (r > 0.85).
- Znotraj zagona je determinističen, seme je treba fiksirati. Celic ne meša, ker fluks računa iz ekspresije celice.

Naredil sem `01_scfea_demo`: 2000 celic P24, preverba orientacije in imen genov, merjenje časa, pregled izhoda.

Odprto: ali so v `data_rna` simboli ali Ensembl ID, imena `module_gene` in `cmMat`. Fluksi scFEA so relativni in odvisni od modela.

## 28. 6. 2026 – Prazne celice v podatkih RNA

Na vzorcu 300 celic P24 (`data_rna_sample.pkl`) je ~30 % praznih celic (kri 34 %, tumor 26 %). Original bere že filtrirane `.mm`, zato v 01 manjka QC.

Po literaturi (Luecken in Theis, 10x, EmptyDrops) sem dodal `min_genes = 200`. Po ponovnem zagonu vseh notebookov so rezultati boljši.

## 4. 7. 2026 – Pregled predobdelave na surovih .h5

- Normalizacija: 01 je delal samo log1p surovih števil (0 do ~11). Avtorji (`data_processing.py:49-58`) naredijo log(x + eps), po celici odštejejo minimum in delijo z maksimumom (izraženi geni ~0.84–1.0, ničle ostanejo 0). To sem implementiral za redko matriko. Razlika z avtorjevo različico je največ 8e-8. Surova števila shranim v `data_rna_counts.pkl` za flukse.
- Ključ (`barcode_clean`, `Patient_ID`) ne loči tkiva in faze: 102 ključa sta v krvi in tumorju, 141 v več fazah. Okoli 283 celic (0.3 %) ima lahko napačno tkivo ali fazo. Nisem popravil.
- V cell-22 sem dodal assert za razpon `data_rna`. Flux demo naj bere `data_rna_counts.pkl`.

## 28. 6. 2026 (zvečer) – TRIM deluje, AUC 0.77

P24 po popravkih (normalizacija, QC): AUC 0.773, osnovni modeli ~0.5 (najboljši knn_rna 0.53). Prej 0.52.

- Spearman med napovedjo in resnico 0.442 (prej 0.019).
- Klonalnost: kri pred rho 0.375, kri po 0.331 (prej 0.079).
- P24 po QC: kri pred 523, kri po 709, tumor pred 1889, tumor po 2796.

## 28. 7. 2026 – AUC 0.87 z avtorjevo normalizacijo in filtrom mito

P24 z avtorjevo normalizacijo in mito % < 20: AUC 0.87 (prej 0.77 z log1p). Najboljši osnovni model knn_both 0.66.

- Youden: recall EXP 0.78 (18/23, prej 0.61), precision 0.51 (prej 0.39), F1 0.62 (prej 0.47). F1 NOT 0.98.
- Spremenil sem dve stvari hkrati. Večji vpliv ima verjetno normalizacija. Filter mito odstrani 2354 celic (1.4 %), mediana mito je 3.7 %.
- 23 pozitivnih, interval zaupanja okoli ±0.1. Za trdnejši sklep rabim več bolnikov.
- ROC je še vedno odrezana zaradi celic z y = 0.
- `pip install --upgrade` je namestil pandas 3, `.pkl` iz pandas 3 se v pandas 2 ne odpre. Brez `--upgrade`, povsod pandas 2.

## 30. 7. 2026 – Prvi izračun fluksov

scFEA ni vzdrževan, zato sem ga prilagodil na pandas 3 (`._append` v `concat`) in vektoriziral pripravo podatkov. Približno 70000 celic na A100: prej 18 h, potem 1 h. Vrednosti so neničelne in na prvi pogled smiselne.

Naprej: kodirnik in dekodirnik za flukse.

## 30. 7. 2026 – Var 1 (fluks + TCR, brez RNA)

Fluks namesto RNA v TRIM (`03_trim_flux`, `04_expansion_eval_flux`).

- Pri prvem treningu je bila izguba rekonstrukcije fluksa 0.000, model je fluks prezrl. Fluksi scFEA so majhni (povprečje ~0.006), izguba TCR pa ~0.15–0.3. RNA tega nima zaradi normalizacije in PCA na visoko variabilnih genih. Rešitev: z-standardizacija po modulih, mean(x²) z 0.0003 na 0.99. Povprečja in std so v `flux_standardize.pkl`.
- P24: fluks + TCR AUC 0.70, RNA + TCR 0.87, NN_both 0.61. Youden: recall EXP 0.70, precision 0.11 (4.4 % ekspandiranih).

Fluks je nad osnovnim modelom in pod RNA + TCR.

## 11. 8. 2026 – Var 2 (RNA + fluks + TCR), pregled kode, priprava ponovitev

P24, en zagon: AUC 0.954, recall EXP 1.0 (23/23), precision 0.40.

Kodo sem pregledal z orodji UI in jo primerjal z originalom. Uhajanja oznak ni.

- Resnica x (velikosti klonov iz `df_all_tcrs`) ne vstopa v izgubo.
- RNA in fluks vplivata na napoved prek skupnega z (povprečje `recon_rna_z`, `recon_flux_z`, `recon_tcr_z`), iz katerega dekodirnik TCR generira psevdoklone.
- Recall 1.0 je posledica velikega deleža y = 0 (filter BB ≥ 2). Hanley-McNeil 95 % IZ [0.894, 1.0], razlika do 0.87 ni značilna (p ~ 0.14). En zagon pri n = 23 in nedeterminističnem treningu.

Ponovitve:

- Seme sem dodal v trening, vrednotenje in flukse (`flux_seed`). Zdaj je ponovljivo.
- Mape `runs/P{pid}/seed{s}/{variant}/`, skupni fluksi po semenih, shranjeni model, izgube in meta.
- Diagnostika za izbor ≥ 5 bolnikov, funkcija za pooled ROC.
- Načrt: 3 variante × ≥ 5 bolnikov × 3 semena.
- Prag psevdoklonov umerjam le na celicah kri pred izpuščenega bolnika, original na `mask_train` vseh bolnikov. Razlog je pomnilnik. Zapisano v cell-25 vseh 03.

## 12. 11. 2026 – Pregled kode in popravki

Pri pisanju diplome sem šel čez vso kodo in jo primerjal z dnevnikom. Podroben seznam je v `diploma/pregled_kode.md`.

Avgust, dopolnjeno iz zapisov zagonov in gita:

- 18. 8.: popravek v rna_tcr (globalni x se je prepisal, napaka v `log_run`).
- 18.–19. 8.: zagoni P24 (rna_tcr semeni 0 in 1, flux_tcr in rna_flux_tcr semena 0–2). rna_tcr P24 seme 2 ni bil pognan.
- 20. 8.: koda za Var 3 (RNA brez genov fluksov + fluks + TCR). Ni bila pognana.
- 21.–22. 8.: zagoni P23, vse tri variante, semena 0–2. flux_tcr seme 0 sem pognal dvakrat (AUC 0.745 in 0.730).
- 31. 8.: zapisi zagonov v `results/logs`.

Skupaj 17 zagonov na dveh bolnikih namesto ≥ 5 bolnikov.

Popravki prejšnjih vnosov:

- 1. 5.: v kodi ni razvrščanja po UMI. Če ima celica več kontigov TRB, ostane prvi v datoteki.
- 1. 5.: `data_rna` je bil log1p le do 4. 7. Končni nabor po QC in brez celic brez TCR ima 69930 celic.
- 13. 5. in 20. 6.: prag se zdaj umeri na TCR, ki jih model generira iz naključnih z s pogoji celic izpuščenega bolnika, ne na vzorcu 100k celic. V 6 od 17 zagonov umerjanje pride do meje 50 ponovitev (prag 0.769), takrat prag ni umerjen.
- 20. 6.: notebook 04 je v gitu od 26. 6.
- 27. 6.: scFEA sem pognal po kosih po 5000 celic (14 kosov). Vsak kos ima svojo mrežo in svojo normalizacijo po celicah, zato je fluks celice odvisen od kosa. Vhod so samo geni modulov.
- 28. 6. (zvečer): vnos sem dodal 28. 7., skupaj z vnosom 4. 7. Filter `min_genes = 200` je v gitu od 4. 7.
- 30. 7.: PCA je na vseh genih, ne na visoko variabilnih. Fluksi so standardizirani na vseh celicah, tudi na izpuščenem bolniku. `flux_standardize.pkl` se ob vsakem zagonu prepiše.
- 30. 7. in 11. 8.: osnovni modeli na fluksih (`knn_flux`, `NN_flux`) so v zapisih zagonov, kode zanje v repozitoriju ni. Shranjeni `04_expansion_eval_flux` ima starejšo različico z `x_rna`.
- 11. 8.: `flux_seed` samo poimenuje izhodno mapo. scFEA ima v kodi `torch.manual_seed(16)`, zato se trije nabori fluksov razlikujejo le zaradi nedeterminizma GPU. Ponovljivo ni: isti zagon je dal 0.745 in 0.730, zagoni z različnimi semeni pa se značilno razlikujejo v 9 od 16 primerjav (DeLong, p < 0.05).
- 11. 8.: AUC 0.954 se ni ponovil. rna_flux_tcr P24 s semeni 0–2 je dal 0.785 do 0.848.
