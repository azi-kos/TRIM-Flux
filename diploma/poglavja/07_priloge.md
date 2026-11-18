<!-- Priloge (neobvezne; ne štejejo v obseg 30 strani ±20 %). Tabela P1 je bila prej tabela 1 v Metodah
     (premaknjena po prvem krogu ocene, 17. 11. 2026). Vir: data_labels.pkl in df_all_tcrs.pkl z Drive (27. 7. 2026). -->

# Priloge

Table: Tabela P1: Število celic po bolnikih in kvadrantih po kontroli kakovosti. GEO je oznaka bolnika v zbirki GEO. Kloni so različna zaporedja CDR3β pri bolniku. Zadnja stolpca navajata celice in klone v krvi pred zdravljenjem, ki so po zdravljenju ekspandirani (razdelek 3.2). Ker smo celice klonov šteli prek vseh bolnikov, imata nekaj ekspandiranih celic tudi bolnika 16 in 20, ki nimata vzorca krvi po zdravljenju. Bolnika, ki smo ju izpustili pri učenju, sta označena krepko.

| Bolnik | GEO | Kri pred | Kri po | Tumor pred | Tumor po | Skupaj | Kloni | Eksp. celice | Eksp. kloni |
|--------|------|-------:|-------:|-------:|-------:|-------:|-------:|--------:|--------:|
| 0 | P01 | 571 | 449 | 0 | 0 | 1020 | 911 | 18 | 7 |
| 1 | P02 | 341 | 369 | 0 | 0 | 710 | 642 | 10 | 4 |
| 2 | P04 | 846 | 800 | 0 | 0 | 1646 | 1395 | 41 | 8 |
| 3 | P05 | 676 | 683 | 0 | 0 | 1359 | 1058 | 75 | 8 |
| 4 | P08 | 447 | 460 | 0 | 0 | 907 | 779 | 26 | 5 |
| 5 | P09 | 854 | 788 | 0 | 0 | 1642 | 1457 | 22 | 12 |
| 6 | P10 | 804 | 473 | 0 | 0 | 1277 | 1118 | 14 | 5 |
| 7 | P12 | 1021 | 1150 | 0 | 0 | 2171 | 1888 | 76 | 22 |
| 8 | P13 | 1088 | 1219 | 0 | 321 | 2628 | 2319 | 2 | 2 |
| 9 | P14 | 409 | 503 | 0 | 1088 | 2000 | 1358 | 54 | 4 |
| 10 | P15 | 842 | 1056 | 0 | 2034 | 3932 | 2583 | 49 | 13 |
| 11 | P16 | 334 | 401 | 0 | 1985 | 2720 | 1553 | 26 | 6 |
| 12 | P17 | 380 | 343 | 0 | 0 | 723 | 709 | 0 | 0 |
| 13 | P18 | 803 | 0 | 2942 | 1781 | 5526 | 2161 | 0 | 0 |
| 14 | P19 | 468 | 421 | 0 | 1602 | 2491 | 2202 | 7 | 6 |
| 15 | P20 | 824 | 537 | 0 | 2425 | 3786 | 2158 | 2 | 2 |
| 16 | P21 | 464 | 0 | 0 | 1950 | 2414 | 1597 | 1 | 1 |
| 17 | P22 | 540 | 494 | 0 | 1915 | 2949 | 2268 | 0 | 0 |
| 18 | P23 | 600 | 918 | 33 | 1600 | 3151 | 2263 | 102 | 24 |
| 19 | P24 | 674 | 467 | 1338 | 464 | 2943 | 2389 | 5 | 4 |
| 20 | P25 | 1196 | 0 | 0 | 604 | 1800 | 1566 | 2 | 2 |
| 21 | P26 | 636 | 228 | 0 | 778 | 1642 | 1376 | 1 | 1 |
| 22 | P27 | 0 | 0 | 1462 | 1379 | 2841 | 1828 | 0 | 0 |
| **23** | **P28** | **1159** | **1290** | **0** | **1811** | **4260** | **2384** | **364** | **27** |
| **24** | **P29** | **523** | **709** | **1888** | **2791** | **5911** | **3947** | **23** | **7** |
| 25 | P30 | 445 | 311 | 0 | 597 | 1353 | 1062 | 8 | 2 |
| 26 | P31 | 580 | 723 | 0 | 1695 | 2998 | 1623 | 133 | 19 |
| 27 | P32 | 438 | 550 | 1039 | 1103 | 3130 | 2194 | 1 | 1 |
| skupaj | | 17 963 | 15 342 | 8702 | 27 923 | 69 930 | | 1062 | 192 |

Table: Tabela P2: Hiperparametri avtokodirnika TCR (razdelek 3.3), scFEA (razdelek 3.4) in modela TRIM (razdelek 3.5). Pri modelu TRIM so enaki za vse različice in bolnike (datoteke `args.txt`). Vsi uporabljajo optimizator Adam.

| Hiperparameter | Vrednost |
|---|---:|
| **avtokodirnik TCR** | |
| število kanalov konvolucijskih plasti (jedro 3, korak 2) | 512, 256, 128 |
| razsežnost vložitve | 100 |
| število korakov učenja | 50 100 |
| velikost paketa | 1024 |
| hitrost učenja | 0,001 |
| utež kazni za velikost vložitve | 0,001 |
| **scFEA** | |
| število epoh (vse celice kosa v vsaki epohi) | 100 |
| hitrost učenja | 0,008 |
| seme v izvorni kodi | 16 |
| največje število celic v kosu | 5000 |
| **model TRIM** | |
| število korakov učenja | 20 000 |
| velikost paketa | 4096 |
| hitrost učenja | 0,001 |
| razsežnost latentnega prostora $d$ | 1024 |
| razsežnost vektorja pogoja | 128 |
| število celic za vektorje pogojev | 5120 |
| utež $\lambda_{\mathrm{KL}}$ | 15 |
| uteži $\lambda_{\mathrm{RNA}}$, $\lambda_{\mathrm{fluks}}$, $\lambda_{\mathrm{TCR}}$ | 1 |
| utež $\lambda_{\mathrm{vl}}$ | 10 |
| rob kontrastne izgube $\delta$ | 10 |

<!-- ZVEZEK: qc_hist_geni.csv in qc_hist_mito.csv → slika P1 (porazdelitev števila izraženih genov na celico z označenim
     pragom 200 in porazdelitev deleža mitohondrijskih štetij z označenim pragom 20 %). -->

**[Slika P1 (porazdelitev števila izraženih genov na celico in deleža mitohondrijskih štetij s pragovoma kontrole kakovosti): vstavim po zagonu statističnega zvezka.]**

<!-- ZVEZEK: umap_bolnik23.csv → slika P2. -->

**[Slika P2 (UMAP resničnih in generiranih celic bolnika 23 v vseh štirih kvadrantih, različica RNA + TCR, seme 0): vstavim po zagonu statističnega zvezka.]**
