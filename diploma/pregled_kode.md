# Pregled kode in primerjava z dnevnikom (12. 11. 2026)

Pregledani so vsi zvezki v `notebooks/`, `notebooks_flux/`, `notebooks_rnaflux/` in `analize/`, ter `trim.py`, `analysis/HNSCC/helpers.py`, `notebooks_flux/experiment_logger.py` in zapisi zagonov v `results/logs/logs/`. Koda ni bila spremenjena. Popravki dnevnika so v vnosu z dne 12. 11. 2026 v `lab_diary.md`.

## A. Vpliva na rezultate ali na opis metod v diplomi

| # | Kaj je v kodi | Kaj piše v dnevniku ali drugje | Posledica za diplomo |
|---|---|---|---|
| 1 | `flux_seed` v `01_scfea_flux` le poimenuje izhodno mapo. scFEA ima v svoji kodi fiksno `torch.manual_seed(16)`. | 11. 8.: seme dodano tudi fluksom, "zdaj reproducibilno". | Trije nabori fluksov se razlikujejo le zaradi nedeterminizma GPU. Isti zagon (flux_tcr, P23, seme 0) je dal AUC 0,745 in 0,730. Že opisano v Metodah (3.5) in Razpravi (5.2, 5.4). |
| 2 | scFEA teče po kosih po 5000 celic (14 kosov). Vsak kos ima svojo mrežo in svojo normalizacijo po celicah. | 27. 6.: "mešanje celic nemogoče". | Fluks celice je odvisen od kosa, v katerega pade. Že v Metodah (3.5) in Razpravi (5.2). Delež variance, ki ga pojasni kos (eta²), izračuna `colab_statistika.ipynb`. |
| 3 | Vhod v scFEA so le geni modulov. | Komentar v zvezku: rezultat "identičen". | Skupni signal celice scFEA izračuna le iz modulnih genov. Ni enako kot pri polnem transkriptomu. Že v Metodah (3.5). |
| 4 | Prag psevdoklonov τ se umeri na TCR, generiranih iz naključnih z s pogoji celic kri pred izpuščenega bolnika (original: vse učne celice). Meja je 50 iteracij. | 13. 5. (vzorec 100k), 11. 8. (le izpuščeni bolnik). | V 6 od 17 zagonov umerjanje doseže mejo (τ = 0,769), prag torej ni umerjen (Razprava 5.2, 5.4). Opis v Metodah (3.7) sem natančneje napisal. |
| 5 | Fluksi so z-standardizirani na vseh celicah, tudi na izpuščenem bolniku. PCA na RNA je izračunan na vseh celicah (kot v originalu). | 30.7.: PCA "na visoko-variabilnih genih". | PCA je na vseh genih. Standardizacija uporabi statistike izpuščenega bolnika (ne oznak). Že v Razpravi (5.4). |
| 6 | Osnovni modeli s fluksi (`knn_flux`, `NN_flux`) so v zapisih zagonov, shranjeni `04_expansion_eval_flux` pa uporablja nedefiniran `x_rna` in ključe rna/tcr/both. | Ni omenjeno. | Koda, ki je dala osnovne modele s fluksi, ni v repozitoriju. V shranjenem zvezku se fluksi berejo nestandardizirani, ni pa znano, kako je bilo v kodi, ki je dala rezultate. V Metodah (3.7) je opomba PREVERI. |
| 7 | Osnovni modeli se testirajo na celicah z novimi TCR (OOD) brez podvojenih, TRIM na vseh celicah kri pred. `MLPClassifier` nima `random_state`. | 20. 6. | Primerjava ni parna (že v Metodah 3.7 in Razpravi 5.4). NN se med zagoni z istim semenom razlikuje. |
| 8 | Oznaka ekspanzije sešteje celice z istim CDR3β čez vse bolnike (kot v originalu). | Ni omenjeno. | Že v Metodah (3.2) in Razpravi (5.4). |
| 9 | Pri več TRB kontigih celice ostane prvi v datoteki (`drop_duplicates` brez razvrščanja). | 1. 5.: razvrščanje po UMI. | Že v Metodah (3.2) in Razpravi (5.4). |
| 10 | Ključ (barkoda, bolnik) ne loči tkiva in faze. | 4. 7.: ~283 celic (0,3 %), nerešeno. | Ostaja. Že v Razpravi (5.4). |
| 11 | Bolnik se celici določi z `groupby('barcode_clean')['Patient_ID'].first()` po vseh vzorcih istih bolnikov. | Ni omenjeno. | Ni preverjeno. Če se ista barkoda pojavi pri dveh bolnikih iz istega nabora, dobi prvega. Preveri se lahko s štetjem takih barkod v metapodatkih. |

## B. Neskladja v dnevniku brez vpliva na rezultate

- Vnosa 4. 7. in 28. 6. (zvečer) sta bila v dnevnik dodana 28. 7. Filter `min_genes = 200`, na katerega se sklicuje vnos 28. 6. (zvečer), je v git od 4. 7.
- Notebook 04 (vnos 20. 6.) je v git od 26. 6.
- Vnos 11. 8. načrtuje ≥ 5 bolnikov in 45 zagonov. Pognanih je 17 zagonov na dveh bolnikih (P23, P24). Manjka rna_tcr P24 seme 2. Var 3 (koda 20. 8.) ni bila pognana.
- AUC 0,954 (11. 8.) se ni ponovil: rna_flux_tcr P24 s semeni 0 do 2 je dal 0,785 do 0,848.
- Končni nabor po QC in izločitvi celic brez TCR ima 69 930 celic (66 829 učnih in 3101 izpuščenih pri P23). V dnevniku ta številka ni zapisana.
- Za zagone 18. do 22. 8., kodo za Var 3 (20. 8.) in shranjevanje zapisov (31. 8.) vnosov ni.

## C. Manjše stvari v kodi (ne vplivajo na rezultate v diplomi)

- Diagnostična celica za QC v zvezkih 04 bere že filtrirane podatke, metrika "RNA kakovost" uporablja naključni z.
- `02_check_flux` ima zastarelo pot.
- `exp_klonov` v diagnostiki šteje celice, ne klonov.
- Markdown v `03_trim` piše, da ima izguba 4 člene (v kodi jih je 5, s fluksi 6).
- Zanka učenja teče 20 001 korak.
- `flux_standardize.pkl` se prepiše ob vsakem zagonu.
- `pip install --upgrade` je še v vseh zvezkih (dnevnik 28. 7. svetuje, naj ga ne bo).
- `analize/01_monte_carlo_analiza.ipynb` združi ocene iz različnih zagonov v en ROC, čeprav ocene med zagoni niso na isti lestvici. Za diplomo se uporablja `diploma/skripte/rezultati_napovedi.py`, ki računa AUC po zagonih.
