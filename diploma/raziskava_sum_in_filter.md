# Šum v RNA in fluksih: ali se splača smer fluksov ali filtra genov (17. 11. 2026)

Vprašanje: ali fluksi, izračunani iz RNA, zmanjšajo šum in povečajo informativnost, ali pa le dodajo napake ocene? Ali je cVAE (TRIM) že dovolj dober in bi bil boljši filter genov iz csGEM (ftINIT) brez fluksov?

## 1. Kakovost informacije: RNA proti fluksom

- **RNA ni izmerjena natančno.** Kapljični scRNA-seq zajame majhen del molekul mRNA v celici, štetja so redka. To posebej velja za šibko izražene gene, med katerimi je veliko encimov. Zaradi tega METAFlux flukse za scRNA-seq računa le na gručah [6 v raziskava_fluksi.md], ftINIT pa za stabilen model potrebuje tisoče združenih celic [5 v raziskava_fluksi.md].
- **Fluksi iz RNA ne morejo vsebovati več informacije o ekspanziji kot RNA, iz katere so izračunani** (neenakost obdelave podatkov, Cover in Thomas, izrek 2.8.1). Lahko pomagajo le kot dober prior (struktura mreže), ki modelu olajša učenje.
- **mRNA slabo napove flukse.** Metode, ki transkriptom vključijo v modele omejitev, niso zanesljivo boljše od pFBA, ki izražanja sploh ne uporabi [1]. Na flukse bolj kot ravni encimov vplivajo koncentracije metabolitov (več kot dvakrat večji vpliv) [2]. Pri celicah T CD4 metabolni preklop ob aktivaciji uravnava translacija mRNA, ki so v mirujočih celicah že prisotne [3].
- **Ocene fluksov dodajo lastne napake:** pravila GPR, pragovi, ciljna funkcija, gojišče. Različne metode se po celicah ne ujemajo (scFEA in Compass r od −0,26 do 0,04) [4].
- **Manj vhodov ni samo po sebi manj šuma.** TRIM RNA že skrči na 100 komponent PCA. Pri 15 modelih z biološkimi potmi so naključno premešane poti dale približno enako dobre rezultate, zato je prednost verjetno posledica redkosti povezav in ne biologije [5]. Pri scRNA-seq je izbira genov pomembna: visoko variabilni geni se obnesejo dobro, naključno izbrani geni pa slabo ohranijo biološke razlike [6].

## 2. Kaj kažejo tvoji zagoni (skripta `diploma/skripte/analiza_psevdoklonov.py`)

Ocena y je 0, če je celica pred zdravljenjem v psevdoklonu z manj kot dvema celicama.

| Bolnik | Različica | AUC zagonov | AUC samo iz y > 0 | AUC znotraj y > 0 |
|---|---|---|---|---|
| 23 | RNA + TCR | 0,906 / 0,782 / 0,923 | 0,828 / 0,844 / 0,827 | 0,772 / 0,224 / 0,821 |
| 24 | RNA + TCR | 0,916 / 0,773 | 0,926 / 0,777 | 0,337 / 0,000 |
| 23 | RNA + fluksi + TCR | 0,742 / 0,909 / 0,909 | 0,816 / 0,835 / 0,837 | 0,254 / 0,783 / 0,764 |
| 24 | RNA + fluksi + TCR | 0,848 / 0,806 / 0,785 | 0,853 / 0,804 / 0,779 | 0,317 / 0,579 / 0,677 |

- Večino AUC prinese že delitev na y = 0 in y > 0, torej ali je celica pred zdravljenjem v psevdoklonu z vsaj dvema celicama. Razvrstitev znotraj y > 0 je nestabilna in včasih slabša od naključja.
- Pri bolniku 23 je med celicami z y = 0 v vseh treh semenih RNA + TCR ekspandiranih le 1 %, med celicami z y > 0 vsaj enkrat pa 53 %.
- Množice celic z y > 0 se med zagoni z RNA in z RNA + fluksi pri bolniku 23 prekrivajo 75–83 % (Jaccard). Fluksi skoraj ne spremenijo, katere celice so v psevdoklonih.
- Delež zagonov, v katerih ima celica y > 0, da AUC 0,927 (bolnik 23) in 0,972 (bolnik 24). Povprečje rangov treh semen RNA + TCR da 0,936 (bolnik 23), povprečje dveh semen pri bolniku 24 pa 0,958. To je bolje od vseh posameznih zagonov RNA + TCR. Številke so raziskovalne, brez intervalov zaupanja.

Razlaga: napovedna moč modela TRIM pri tej oznaki pride predvsem iz klonalnosti, ki jo model razbere iz zaporedij TCR. Podrobnosti transkriptoma (in s tem fluksov) imajo manjšo vlogo. Oznaka ekspanzije primerja surova števila celic klona po zdravljenju in pred njim, ne deležev. Klon z eno celico pred zdravljenjem je ekspandiran le, če ima po zdravljenju vsaj dve celici, kar je pri redkih klonih malo verjetno. Oznaka je zato verjetno močno povezana z velikostjo klona pred zdravljenjem in s številom celic, zajetih ob vsakem času (bolnik 24: 523 celic pred, 709 po zdravljenju). To je treba preveriti z izhodiščem iz prave velikosti klona pred zdravljenjem.

## 3. Priporočilo

1. **Najprej preverba v Colabu (okoli 30 min):** AUC prave velikosti klona v krvi pred zdravljenjem (iz `df_all_tcrs` in `data_labels`) za bolnika 23 in 24 ter AUC modela TRIM znotraj skupin klonov enake velikosti. Če izhodišče doseže ~0,9, vhodni podatki (RNA, fluksi, filter) na AUC skoraj ne vplivajo.
2. **Filter genov namesto novih fluksov kot prvi poskus:** csGEM iz ftINIT (vse celice T, brez bolnikov 23 in 24), nato TRIM z RNA samo iz genov csGEM in z RNA brez teh genov (koda za Var 3 že obstaja). To neposredno pove, ali metabolni geni nosijo signal, brez napak ocene fluksov, in je zgornja referenca za vsak fluks iz istih genov.
3. **Novi fluksi (kvadratni program v csGEM) le, če metabolni geni dodajo signal** onkraj velikosti klona.
4. **Pri vseh različicah združevati semena** (povprečje rangov), ker je razlika med semeni večja od razlike med različicami.

## Viri

1. Machado D, Herrgård M. Systematic evaluation of methods for integration of transcriptomic data into constraint-based models of metabolism. PLoS Comput Biol 2014;10(4):e1003580. https://doi.org/10.1371/journal.pcbi.1003580
2. Hackett SR in sod. Systems-level analysis of mechanisms regulating yeast metabolic flux. Science 2016;354(6311):aaf2786. https://doi.org/10.1126/science.aaf2786
3. Ricciardi S in sod. The translational machinery of human CD4+ T cells is poised for activation and controls the switch from quiescence to metabolic remodeling. Cell Metab 2018;28(6):895–906.e5. https://doi.org/10.1016/j.cmet.2018.08.009
4. Amiss JMB, Jabbari H, Lum J. A method to identify high consensus predictions of single-cell metabolic flux. bioRxiv 2024. https://doi.org/10.1101/2024.01.15.572211
5. Caranzano I in sod. Sparsity is all you need: rethinking biological pathway-informed approaches in deep learning. arXiv 2025:2505.04300. https://arxiv.org/abs/2505.04300
6. Zappia L in sod. Feature selection methods affect the performance of scRNA-seq data integration and querying. Nat Methods 2025;22:834–844. https://doi.org/10.1038/s41592-025-02624-3
7. Cover TM, Thomas JA. Elements of Information Theory. 2. izd. Wiley, 2006 (izrek 2.8.1, neenakost obdelave podatkov).
