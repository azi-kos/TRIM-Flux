# Preverba: koliko AUC modela TRIM razloži velikost klona (17. 11. 2026)

Vprašanje: ali visok AUC modela TRIM pride iz transkriptoma (in fluksov) ali že iz velikosti klona pred zdravljenjem, ki jo model razbere iz zaporedij TCR?

Zvezek `diploma/skripte/colab_preverba_klonov.ipynb`, pognan na `data_labels.pkl` in `df_all_tcrs.pkl` z Drive (različica z dne 27. 7. 2026) in na 17 zapisih zagonov iz `results/logs` (enaki kot na Drive). Oznake se ujemajo z `x_true` v vseh zapisih. Rezultati so v `diploma/rezultati/preverba_klonov/`. Intervali so 95-% intervali bootstrap po klonih (2000 ponovitev), ker celice istega klona niso neodvisne.

## 1. Velikost klona sama doseže AUC modela TRIM

`krv_pred` je število celic z istim CDR3β v krvi pred zdravljenjem. Ansambel je povprečje rangov ocen vseh semen iste različice.

| | Bolnik 23 | Bolnik 24 |
|---|---|---|
| Celice (ekspandirane) | 1159 (364) | 523 (23) |
| Kloni (ekspandirani) | 634 (27) | 497 (7) |
| Velikost klona (`krv_pred`) | 0,958 (0,746–0,984) | 0,950 (0,721–0,998) |
| RNA + TCR, ansambel | 0,936 (0,791–0,969) | 0,958 (0,832–0,992) |
| RNA + fluksi + TCR, ansambel | 0,929 (0,784–0,964) | 0,872 (0,446–0,972) |
| Fluksi + TCR, ansambel | 0,796 (0,545–0,866) | 0,916 (0,736–0,971) |
| Vsi zagoni, ansambel | 0,940 (0,769–0,975) | 0,977 (0,913–0,995) |
| Posamezni zagoni TRIM | 0,476–0,923 | 0,659–0,916 |

- Noben posamezen zagon ni boljši od velikosti klona. Ansambli so enako dobri: razlika vseh zagonov proti `krv_pred` je −0,018 (−0,048 do 0,028) pri bolniku 23 in 0,027 (−0,013 do 0,193) pri bolniku 24.
- Ekspandirane celice izvirajo iz malo klonov. Pri bolniku 23 je 309 od 364 ekspandiranih celic iz štirih klonov z več kot 10 celicami, pri bolniku 24 je 23 celic iz sedmih klonov. AUC na celicah zato v resnici sloni na 27 oziroma 7 klonih.
- Med kloni z vsaj dvema celicama je velikost klona boljša od modela TRIM (bolnik 23: 0,925 proti 0,66–0,86 za ansamble).
- Če celice razvrstimo po velikosti klona in vezi razbijemo z oceno TRIM, se AUC pri različicah z RNA nekoliko poveča (bolnik 23 z 0,958 na 0,965–0,966, interval razlike 0,002–0,05).

## 2. Kje TRIM doda nekaj: kloni z eno celico pred zdravljenjem

Pri bolniku 23 je 548 klonov z eno celico v krvi pred zdravljenjem, ekspandiranih (po zdravljenju vsaj dve celici) je 10. Velikost klona jih ne loči (AUC 0,5).

| Bolnik 23, kloni z eno celico | AUC |
|---|---|
| RNA + fluksi + TCR, ansambel | 0,917 (0,866–0,961) |
| RNA + TCR, ansambel | 0,849 (0,701–0,969) |
| Fluksi + TCR, ansambel | 0,443 (0,226–0,652) |

Na ravni klonov (vsak klon šteje enkrat, parni DeLong) je pri bolniku 23 velikost klona 0,766, RNA + TCR 0,890 (p = 0,005), RNA + fluksi + TCR 0,895 (p = 0,005), fluksi + TCR pa 0,558 (p = 0,001, slabše). Bolnik 24 ima le dva ekspandirana klona z eno celico, zato tam tega ni mogoče preveriti (ravni klonov: velikost klona 0,846, RNA + TCR 0,966, p = 0,20).

Signal onkraj velikosti klona je torej samo v različicah z RNA. Fluksi brez RNA ga nimajo, fluksi ob RNA ga na desetih klonih ne povečajo značilno.

## 3. Oznaka je odvisna od globine vzorčenja

- Med 24 bolniki z vzorci krvi pred zdravljenjem in po njem je delež ekspandiranih celic povezan z razmerjem števila zajetih celic po/pred (Spearman 0,58, p = 0,003). Pri oznaki iz deležev (delež klona med celicami bolnika ob tem času) je povezava šibkejša (0,32, p = 0,13).
- Pri bolniku 23 oznaka iz deležev spremeni 62 od 1159 celic (364 → 302 ekspandiranih), AUC se malo zniža (velikost klona 0,931, RNA + TCR 0,894). Pri bolniku 24 se spremeni ena celica.
- Združevanje števil čez bolnike (A9 v `pregled_kode.md`) na izpuščena bolnika skoraj ne vpliva (0 oziroma 1 celica).
- Bolnik 23 je izjemno klonalen: le 47 % celic krvi pred zdravljenjem je v klonih z eno celico (drugi bolniki 63–99 %), ekspandiranih je 31 % celic (drugi 0–23 %).
- Bolnik 23 nima vzorca tumorja pred zdravljenjem, pri bolniku 24 velikost klona v tumorju pred zdravljenjem ekspanzije ne napove (0,535).

## 4. Posledice za nadaljnje delo

1. Skupni AUC na celicah ni primerno merilo za primerjavo vhodov (RNA, fluksi, filter genov), ker ga velikost klona že skoraj zapolni. Merilo za primerjave: AUC znotraj klonov z eno celico in AUC na ravni klonov, vedno ob izhodišču velikosti klona in z intervali po klonih.
2. Moč je majhna: pri izpuščenih bolnikih je 12 ekspandiranih klonov z eno celico (10 + 2). Bolniki z največ ekspandiranimi kloni so 18 (24 klonov, 9 z eno celico), 7 (22, 8), 26 (19, 5) in 10 (13, 6). Čez vseh 24 bolnikov je ekspandiranih klonov z eno celico okoli 70.
3. Filter genov (csGEM) lahko najprej preverimo poceni: preprost klasifikator na celicah krvi pred zdravljenjem iz klonov z eno celico, z izpuščanjem po enega bolnika (vseh 24), z značilkami iz vseh genov, iz metabolnih genov, iz ostalih genov in iz fluksov scFEA. To teče na CPU v nekaj minutah. Nove zagone TRIM (GPU) naredimo le, če metabolni geni pokažejo signal.
4. Za diplomo: v Rezultate tabelo z izhodiščem velikosti klona ob AUC modela TRIM, AUC znotraj klonov z eno celico in na ravni klonov. V Razpravo: AUC na celicah je pretežno posledica velikosti klona, oznaka iz surovih števil je odvisna od globine vzorčenja, ekspandiranih klonov je malo.
