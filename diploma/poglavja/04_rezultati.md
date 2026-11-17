<!-- ~5 strani. Samo rezultati (razprava je v ločenem poglavju, navodila Biokemija v1.2).
     Vse številke so iz zapisov zagonov results/logs/logs/*.json, izračunane s skripto
     diploma/skripte/rezultati_napovedi.py (izhodi v diploma/rezultati/). Številke z oznako »ZVEZEK«
     pridejo iz diploma/skripte/colab_statistika.ipynb.
     Kontrolni seznam: vsak prikazan rezultat opisan v besedilu s sklicem (»slika 6«, »tabela 4«, mala začetnica);
     slike lastne, vse oznake v slovenščini, decimalna vejica; naslov in opis POD sliko, pri tabeli NAD tabelo. -->

# Rezultati

## Kontrola kakovosti in opis podatkov

<!-- ZVEZEK: qc_geni_na_celico.csv in qc_po_datotekah.csv → slika 5 (porazdelitev števila izraženih genov na celico
     z označenim pragom 200 in porazdelitev deleža mitohondrijskih štetij z označenim pragom 20 %), število celic,
     ki jih odstrani vsak korak; tabela1_celice_po_bolnikih.csv → stavek o številu celic in klonov po kvadrantih. -->

**[Besedilo in slika 5 (porazdelitev števila izraženih genov na celico in deleža mitohondrijskih štetij s pragovoma kontrole kakovosti): vpišem po zagonu statističnega zvezka.]**

## Napoved klonske ekspanzije pri izpuščenih bolnikih

Model TRIM smo ovrednotili pri dveh izpuščenih bolnikih, s tremi različicami in tremi semeni, skupaj v 17 zagonih (tabela 4). Pri bolniku 23 smo napoved vrednotili na 1159 celicah iz krvi pred zdravljenjem, od katerih jih je bilo 364 iz ekspandiranih klonov. Pri bolniku 24 smo jo vrednotili na 523 celicah, od katerih jih je bilo 23 iz ekspandiranih klonov. Krivulje ROC vseh zagonov prikazuje slika 6, AUC z intervali zaupanja in AUC izhodiščnih klasifikatorjev pa slika 7.

Table: Tabela 4: AUC napovedi klonske ekspanzije modela TRIM pri izpuščenih bolnikih 23 in 24. V oklepaju je 95-% interval zaupanja po DeLongu; »–« pomeni, da zagona nismo izvedli. Povprečje je izračunano po razpoložljivih semenih.

| Različica | Seme | Bolnik 23 | Bolnik 24 |
|----------------------|-------------|--------------------------|--------------------------|
| RNA + TCR | 0 | 0,906 (0,887–0,925) | 0,916 (0,856–0,976) |
|  | 1 | 0,782 (0,754–0,809) | 0,773 (0,671–0,876) |
|  | 2 | 0,923 (0,905–0,941) | – |
|  | povprečje | 0,870 | 0,845 |
| fluksi + TCR | 0 | 0,745 (0,715–0,774) | 0,813 (0,720–0,906) |
|  | 1 | 0,476 (0,442–0,510) | 0,659 (0,567–0,750) |
|  | 2 | 0,780 (0,752–0,808) | 0,851 (0,779–0,922) |
|  | povprečje | 0,667 | 0,774 |
| RNA + fluksi + TCR | 0 | 0,742 (0,713–0,771) | 0,848 (0,757–0,939) |
|  | 1 | 0,908 (0,889–0,928) | 0,806 (0,706–0,907) |
|  | 2 | 0,909 (0,890–0,928) | 0,785 (0,681–0,889) |
|  | povprečje | 0,853 | 0,813 |

Model je velikemu deležu celic pripisal oceno $y_i = 0$ (enačba 7), ker je bila povprečna velikost njihovega psevdoklona v krvi pred zdravljenjem manjša od 2. Pri različicah z RNA je bil ta delež pri bolniku 23 med 48 in 54 %, pri bolniku 24 pa med 90 in 96 %. Celic z oceno 0 model med seboj ne razvrsti, zato imajo krivulje ROC pri bolniku 24 dolg raven odsek do točke (1, 1) (slika 6B).

![Slika 6: Krivulje ROC napovedi klonske ekspanzije pri izpuščenem bolniku 23 (A) in 24 (B). Vsaka črta je en zagon (seme); barva in vrsta črte označujeta različico modela. Diagonala ustreza naključnemu ugibanju.](slike/roc_izpusceni_bolniki.png){width=100%}

![Slika 7: AUC napovedi klonske ekspanzije po različicah in semenih pri izpuščenem bolniku 23 (A) in 24 (B). Točka je AUC posameznega zagona, navpična črta njen 95-% interval zaupanja po DeLongu, številka pod točko pa seme. Vodoravna siva črta ob točki je AUC najboljšega izhodiščnega klasifikatorja iz istega zagona (največja vrednost med klasifikatorji kNN in NN na vseh vhodih, tabela 5).](slike/auc_semena_izhodisca.png){width=100%}

Table: Tabela 5: AUC izhodiščnih klasifikatorjev. Pri vhodih, ki smo jih uporabili v več zagonih, je naveden razpon vrednosti.

| Klasifikator | Vhod | Bolnik 23 | Bolnik 24 |
|---|---|---|---|
| kNN | RNA | 0,485 | 0,652 |
| kNN | TCR | 0,518 | 0,419 |
| kNN | RNA in TCR | 0,488 | 0,659 |
| kNN | fluksi | 0,506–0,515 | 0,461–0,486 |
| kNN | fluksi in TCR | 0,464 | 0,206 |
| NN | RNA | 0,472 | 0,348 |
| NN | TCR | 0,476–0,481 | 0,457–0,477 |
| NN | RNA in TCR | 0,397–0,444 | 0,153–0,191 |
| NN | fluksi | 0,541–0,548 | 0,437–0,496 |
| NN | fluksi in TCR | 0,628–0,733 | 0,206–0,503 |

## Ponovitev modela TRIM na transkriptomu in TCR

Različica RNA + TCR je pri bolniku 23 dosegla AUC 0,906, 0,782 in 0,923 (povprečje 0,870), pri bolniku 24 pa 0,916 in 0,773 (povprečje 0,845; tabela 4). V vseh petih zagonih je bila spodnja meja 95-% intervala zaupanja nad 0,5 in nad AUC najboljšega izhodiščnega klasifikatorja (slika 7). Izhodiščni klasifikatorji na RNA, TCR ali obojem so pri bolniku 23 dosegli AUC med 0,397 in 0,518, pri bolniku 24 pa med 0,153 in 0,659 (tabela 5). Najboljši izhodiščni klasifikator je bil pri bolniku 23 kNN na TCR (AUC 0,518), pri bolniku 24 pa kNN na RNA in TCR (AUC 0,659). Pri bolniku 24 je bila spodnja meja intervala zaupanja modela TRIM s semenom 1 (0,671) le malo nad to vrednostjo. Večplastni perceptron je imel pri bolniku 24 na RNA ter na RNA in TCR AUC pod 0,5 (0,153–0,348).

Pri pragu z največjim Youdenovim indeksom je različica RNA + TCR pri bolniku 23 kot ekspandirane prepoznala 84–92 % celic iz ekspandiranih klonov pri specifičnosti 0,78–0,94 in preciznosti 0,66–0,87. Pri bolniku 24 je prepoznala 91 % (seme 0) oziroma 57 % (seme 1) celic iz ekspandiranih klonov pri specifičnosti 0,94 oziroma 0,99 in preciznosti 0,40 oziroma 0,68.

<!-- ZVEZEK: umap_bolnik23.csv → slika 8. -->

**[Slika 8 (UMAP resničnih in generiranih celic bolnika 23 v vseh štirih kvadrantih, različica RNA + TCR, seme 0): vstavim po zagonu statističnega zvezka.]**

## Metabolni fluksi kot samostojna modaliteta

Različica fluksi + TCR je pri bolniku 23 dosegla AUC med 0,476 in 0,780 (povprečje 0,667), pri bolniku 24 pa med 0,659 in 0,851 (povprečje 0,774; tabela 4). Pri bolniku 23 je bila AUC zagona s semenom 1 (0,476; 95-% interval zaupanja 0,442–0,510) na ravni naključnega ugibanja. Pri bolniku 24 je bila v vseh treh zagonih spodnja meja intervala zaupanja nad 0,5 in nad AUC vseh izhodiščnih klasifikatorjev (največ 0,503).

Izhodiščni klasifikatorji na samih fluksih so dosegli AUC med 0,437 in 0,548 (tabela 5). Večplastni perceptron na fluksih in TCR je pri bolniku 23 dosegel AUC med 0,628 in 0,733, v ponovljenem zagonu s semenom 0 pa 0,830. Spodnja meja intervala zaupanja modela TRIM je bila nad AUC tega klasifikatorja le v zagonu s semenom 2 (slika 7A).

<!-- ZVEZEK: fluksi_ponovljivost.csv → ujemanje fluksov med semeni scFEA (korelacija, največja razlika);
     fluksi_kosi_eta2.csv → delež variance fluksov, ki ga pojasni pripadnost kosu 5000 celic (mediana, razpon po modulih);
     fluksi_moduli_ekspanzija.csv → število modulov s q < 0,05 pri vseh bolnikih ter pri bolnikih 23 in 24,
     največja in najmanjša AUC modula, imena modulov → slika 9. -->

**[Besedilo in slika 9 (ponovljivost fluksov med zagoni scFEA, vpliv kosov celic in razlike v fluksih posameznih modulov med celicami iz ekspandiranih in neekspandiranih klonov): vpišem po zagonu statističnega zvezka.]**

## Združevanje transkriptoma, fluksov in TCR

Različica RNA + fluksi + TCR je dosegla povprečno AUC 0,853 pri bolniku 23 in 0,813 pri bolniku 24, različica RNA + TCR pa 0,870 in 0,845 (tabela 4). Na istih celicah in z istim semenom je bila razlika v AUC med različicama RNA + fluksi + TCR in RNA + TCR v treh od petih parov negativna, v dveh pa pozitivna (tabela 6). Statistično značilni sta bili dve razliki pri bolniku 23, ki sta imeli nasprotna predznaka: s semenom 0 je bila AUC različice s fluksi nižja za 0,164, s semenom 1 pa višja za 0,127 (obe *p* < 0,001). Pri bolniku 24 nobena razlika ni bila statistično značilna (*p* ≥ 0,18). Različica fluksi + TCR je imela v vseh petih parih nižjo AUC kot različica RNA + TCR, z razliko od 0,104 do 0,306 (*p* ≤ 0,020).

Table: Tabela 6: Razlika v AUC (ΔAUC) različic RNA + fluksi + TCR in fluksi + TCR do različice RNA + TCR na istih celicah pri istem bolniku in semenu, s parnim DeLongovim testom. Pozitivna razlika pomeni višjo AUC različice s fluksi; pri bolniku 24 s semenom 2 zagona RNA + TCR nismo izvedli.

| Bolnik | Seme | RNA + fluksi + TCR: ΔAUC | *p* | fluksi + TCR: ΔAUC | *p* |
|----------|--------|------------------------|------------|------------------------|------------|
| 23 | 0 | −0,164 | < 0,001 | −0,161 | < 0,001 |
| 23 | 1 | 0,127 | < 0,001 | −0,306 | < 0,001 |
| 23 | 2 | −0,014 | 0,17 | −0,143 | < 0,001 |
| 24 | 0 | −0,068 | 0,18 | −0,104 | 0,014 |
| 24 | 1 | 0,033 | 0,45 | −0,114 | 0,020 |

## Ponovljivost med semeni in umerjanje praga psevdoklonov

AUC iste različice pri istem bolniku se je med semeni razlikovala do 0,143 pri različici RNA + TCR, do 0,167 pri različici RNA + fluksi + TCR in do 0,304 pri različici fluksi + TCR (tabela 4). Parni DeLongov test je pokazal statistično značilno razliko (*p* < 0,05) med zagoni z različnimi semeni v 9 od 16 parov. Spearmanov koeficient korelacije med ocenami $y_i$ dveh zagonov z različnima semenoma je bil pri različici RNA + TCR med 0,50 in 0,77, pri različici RNA + fluksi + TCR med 0,36 in 0,76, pri različici fluksi + TCR pa med −0,03 in 0,30. Ponovljeni zagon različice fluksi + TCR pri bolniku 23 s semenom 0 je dal AUC 0,730 namesto 0,745 (razlika 0,015; *p* = 0,004), Spearmanov koeficient korelacije med ocenama obeh zagonov pa je bil 0,97.

Umerjanje praga psevdoklonov se je v 6 od 17 zagonov ustavilo po 50 ponovitvah pri vrednosti $\tau = 10 \cdot 0{,}95^{50} \approx 0{,}77$, ne da bi število psevdoklonov doseglo število klonov. To se je zgodilo v vseh treh zagonih različice fluksi + TCR pri bolniku 24, v obeh zagonih različice RNA + TCR pri bolniku 24 in v zagonu različice fluksi + TCR s semenom 1 pri bolniku 23. V ostalih zagonih je bil umerjeni prag pri bolniku 23 med 0,99 in 12,16, pri bolniku 24 pa med 0,85 in 4,40. V zagonu različice fluksi + TCR s semenom 1 pri bolniku 23 je imela vsaka celica oceno $y_i > 0$, AUC tega zagona pa je bila 0,476.

## Velikost klona pred zdravljenjem kot izhodišče

Velikost klona v krvi pred zdravljenjem je ekspanzijo napovedala z AUC 0,958 pri bolniku 23 in 0,950 pri bolniku 24 (tabela 7). To je več kot v katerem koli posameznem zagonu modela TRIM (od 0,476 do 0,923 pri bolniku 23 in od 0,659 do 0,916 pri bolniku 24; tabela 4). Ansambel semen je imel v vseh šestih kombinacijah bolnika in različice višjo AUC od najboljšega posameznega zagona iste različice, pri različici RNA + TCR na primer 0,936 pri bolniku 23 in 0,958 pri bolniku 24. Razlika med ansamblom vseh zagonov in velikostjo klona je bila pri bolniku 23 −0,018 (95-% interval bootstrap od −0,048 do 0,028), pri bolniku 24 pa 0,027 (od −0,013 do 0,193). Število celic klona samo pri izpuščenem bolniku je dalo skoraj enako AUC kot število prek vseh bolnikov (0,959 in 0,929), število celic klona v tumorju pred zdravljenjem pa ekspanzije pri bolniku 24 ni napovedalo (0,535). Bolnik 23 vzorca tumorja pred zdravljenjem ni imel.

Table: Tabela 7: AUC velikosti klona v krvi pred zdravljenjem in ansamblov modela TRIM pri izpuščenih bolnikih. Ansambel je povprečje rangov ocen vseh zagonov različice, »vsi zagoni« pa povprečje rangov vseh zagonov bolnika. Pri celicah je v oklepaju 95-% interval bootstrap po klonih, pri klonih pa 95-% interval zaupanja po DeLongu. Pri bolniku 24 sta bila med kloni z eno celico ekspandirana le dva, zato AUC v tej skupini ne navajamo.

| Bolnik | Ocena | Vse celice | Celice iz klonov z eno celico | Kloni |
|--------|----------------------|---------------------|---------------------|---------------------|
| 23 | število (ekspandiranih) | 1159 (364) | 548 (10) | 634 (27) |
|  | velikost klona | 0,958 (0,746–0,984) | 0,500 | 0,766 (0,670–0,863) |
|  | RNA + TCR | 0,936 (0,791–0,969) | 0,849 (0,701–0,969) | 0,890 (0,829–0,951) |
|  | fluksi + TCR | 0,796 (0,545–0,866) | 0,443 (0,226–0,652) | 0,558 (0,433–0,683) |
|  | RNA + fluksi + TCR | 0,929 (0,784–0,964) | 0,917 (0,866–0,961) | 0,895 (0,858–0,932) |
|  | vsi zagoni | 0,940 (0,769–0,975) | 0,814 (0,609–0,975) | 0,864 (0,785–0,944) |
| 24 | število (ekspandiranih) | 523 (23) | 475 (2) | 497 (7) |
|  | velikost klona | 0,950 (0,721–0,998) | – | 0,846 (0,662–1,000) |
|  | RNA + TCR | 0,958 (0,832–0,992) | – | 0,966 (0,946–0,987) |
|  | fluksi + TCR | 0,916 (0,736–0,971) | – | 0,838 (0,707–0,969) |
|  | RNA + fluksi + TCR | 0,872 (0,446–0,972) | – | 0,683 (0,468–0,898) |
|  | vsi zagoni | 0,977 (0,913–0,995) | – | 0,952 (0,903–1,000) |

Intervali bootstrap po klonih so bili precej širši od DeLongovih intervalov na celicah (pri velikosti klona pri bolniku 23 od 0,746 do 0,984 namesto od 0,945 do 0,972), ker so ekspandirane celice izvirale iz malo klonov. Pri bolniku 23 je bilo 364 celic iz ekspandiranih klonov iz 27 od 634 klonov, od tega 309 celic iz štirih klonov z več kot 10 celicami v krvi pred zdravljenjem. Pri bolniku 24 so 23 celic iz ekspandiranih klonov izvirale iz 7 od 497 klonov. Delež ekspandiranih klonov je naraščal z velikostjo klona. Med kloni z eno celico je bil 1,8 % pri bolniku 23 in 0,4 % pri bolniku 24, vsi kloni z več kot 10 celicami pa so bili ekspandirani (slika 10A, B).

Med celicami iz klonov z eno celico v krvi pred zdravljenjem velikost klona celic ne loči. Pri bolniku 23 je bilo takih klonov 548, ekspandiranih (z vsaj dvema celicama po zdravljenju) pa 10. Ansambel različice RNA + fluksi + TCR je med njimi dosegel AUC 0,917, ansambel različice RNA + TCR 0,849, ansambel različice fluksi + TCR pa 0,443 (tabela 7). Na ravni klonov sta imeli različici z RNA pri bolniku 23 višjo AUC kot velikost klona (0,890 in 0,895 proti 0,766; parni DeLongov test, obakrat *p* = 0,005), različica fluksi + TCR pa nižjo (0,558; *p* = 0,001). Pri bolniku 24 se na ravni klonov nobena različica od velikosti klona ni statistično značilno razlikovala (*p* ≥ 0,085). Med celicami iz klonov z vsaj dvema celicama je bila velikost klona pri bolniku 23 boljša od vseh ansamblov (AUC 0,925 proti 0,661–0,856). Ko smo celice razvrstili po velikosti klona in celice iz klonov enake velikosti še po oceni modela, se je AUC pri bolniku 23 povečala z 0,958 na 0,965 z ansamblom različice RNA + TCR (interval razlike od 0,002 do 0,046) in na 0,966 z ansamblom različice RNA + fluksi + TCR (od 0,002 do 0,051).

![Slika 10: Ekspanzija klonov glede na velikost klona in globino vzorčenja. A, B: delež ekspandiranih klonov med kloni izpuščenega bolnika 23 (A) in 24 (B), razvrščenimi po številu celic v krvi pred zdravljenjem; nad stolpcem je število ekspandiranih in vseh klonov, »–« pomeni, da v skupini ni bilo klonov. C: delež celic iz ekspandiranih klonov v krvi pred zdravljenjem in razmerje med številom celic v krvi po zdravljenju in pred njim pri 24 bolnikih z obema vzorcema krvi; ρ je Spearmanov koeficient korelacije, modri točki sta izpuščena bolnika.](slike/velikost_klona.png){width=100%}

Delež celic iz ekspandiranih klonov se je med 24 bolniki z obema vzorcema krvi gibal od 0 do 31 % in je bil povezan z razmerjem med številom celic, zajetih v krvi po zdravljenju in pred njim (Spearmanov ρ = 0,58, *p* = 0,003; slika 10C). Pri oznaki iz deležev je bila povezava šibkejša (ρ = 0,32, *p* = 0,13). Pri bolniku 23 je oznaka iz deležev spremenila oznako 62 od 1159 celic (302 ekspandirane namesto 364), velikost klona je pri njej dosegla AUC 0,931, ansambel različice RNA + TCR pa 0,894. Pri bolniku 24 se je spremenila oznaka ene celice. Bolnik 23 je bil med vsemi bolniki najbolj klonalen: v klonih z eno celico je bilo le 47 % njegovih celic iz krvi pred zdravljenjem, pri ostalih bolnikih z obema vzorcema krvi pa od 63 do 99 %.
