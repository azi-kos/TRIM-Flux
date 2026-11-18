<!-- Materiali in metode (Navodila 2026 + Biokemija v1.2): dovolj natančno za ponovitev, logični vrstni red,
     vsaka metoda iz literature z virom, verzije programske opreme in strojna oprema.
     Vse številke so iz kode (notebooks*/), results/*/args.txt, results/*/heldout_meta.json, results/logs/logs/*.json
     in lab_diary.md. Številke z oznako »ZVEZEK« pridejo iz diploma/skripte/colab_statistika.ipynb.
     Zgradba po prvem krogu ocene (17. 11. 2026): modalitete (3.2–3.4) pred modelom (3.5), napoved (3.6),
     vrednotenje (3.7), statistika (3.8), naknadna analiza velikosti klona (3.9). Tabela celic po bolnikih je v Prilogi (P1), hiperparametri tudi (P2). -->

# Materiali in metode

Potek dela od surovih podatkov do vrednotenja prikazuje slika 3. Vso analizo smo izvedli v zvezkih Jupyter v okolju Google Colab; koda je v repozitoriju TRIM-Flux (razdelek 3.10).

![Slika 3: Potek dela od surovih podatkov do vrednotenja napovedi. Iz istih celic smo pripravili tri modalitete (transkriptom, TCR in metabolne flukse), ki smo jih v različnih kombinacijah dali v model TRIM. Številke v oklepajih so razdelki in enačbe tega poglavja.](slike/potek_dela.png){width=100%}

## Podatki

Uporabili smo javno dostopne podatke raziskave Luome in sodelavcev [@luoma2022] iz zbirke Gene Expression Omnibus (GEO) z oznako GSE200996 [@barrett2013]. Raziskava je zajela bolnike s HNSCC, ki so pred operacijo prejeli neoadjuvantno imunoterapijo [@luoma2022]. Iz zbirke smo prenesli matrike štetij prepisov (datoteke `.h5`), anotacije kontigov TCR (datoteke `filtered_contig_annotations`) in metapodatke celic T CD4+ in CD8+ iz krvi in tumorja. Metapodatki za vsako celico navajajo bolnika, tkivo, podtip in oznako vzorca.

Vzorce z oznakama B1 (kri) in Pre-Tx (tumor) smo obravnavali kot vzorce pred zdravljenjem, vzorce z oznakama B2 (kri) in Post-Tx (tumor) pa kot vzorce po zdravljenju. Vzorcev z oznako B3 nismo uporabili, ker zanje zbirka ne vsebuje matrik štetij. Kombinacijo tkiva (kri ali tumor) in časa vzorčenja (pred zdravljenjem ali po njem) imenujemo kvadrant. Bolnike smo oštevilčili od 0 do 27 po abecednem vrstnem redu njihovih oznak v zbirki GEO. V nadaljevanju jih navajamo s temi indeksi, oznake GEO in število celic po kvadrantih pa so v tabeli P1 v Prilogi.

## Predobdelava podatkov ter določitev klonov in ekspanzije

Podatke smo predobdelali s paketom scanpy [@wolf2018] po splošnih priporočilih za analizo podatkov scRNA-seq [@luecken2019]. Datoteke `.h5` vsebujejo surove matrike štetij, v katerih so poleg celic tudi prazne mikrokapljice in poškodovane celice. Kontrolo kakovosti smo zato izvedli za vsako datoteko posebej v treh korakih. Najprej smo obdržali celice z vsaj 200 izraženimi geni; prag smo izbrali na podlagi porazdelitve števila izraženih genov na celico (slika P1 v Prilogi). Nato smo obdržali samo celice, ki so v metapodatkih avtorjev [@luoma2022]. Avtorji so vzorce več bolnikov sekvencirali skupaj, zato se ista črtna koda lahko pojavi pri več bolnikih, in celico smo vedno označili s parom (črtna koda, bolnik). Na koncu smo odstranili celice z 20 % ali več štetij mitohondrijskih genov (imena z začetkom MT-), ker visok delež kaže na poškodovano celico [@luecken2019].

Iz anotacij kontigov TCR smo obdržali verige β (TRB), ki so bile produktivne, visoko zanesljive in so imele določeno zaporedje CDR3. Če je imela celica več takšnih kontigov, smo obdržali prvega. Izločili smo celice brez zaporedja CDR3β. Po vseh korakih je ostalo 69 930 celic s 36 601 genom (razdelek 4.1).

<!-- ZVEZEK: qc_po_datotekah.csv za število celic, ki jih odstrani vsak korak. -->

Par (črtna koda, bolnik) ne loči tkiva in časa vzorčenja istega bolnika: 102 para v metapodatkih kažeta na kri in tumor hkrati, 141 pa na več časov vzorčenja. Obdržali smo prvo oznako, zato ima lahko napačno oznako tkiva ali časa največ okoli 0,3 % celic.

Štetja smo normalizirali enako kot avtorji modela TRIM [@he2026trim]. Za vsako celico $i$ smo štetja $x_{ig}$ logaritmirali in jih preslikali na interval od 0 do 1:

$$\tilde{x}_{ig} = \frac{\ln(x_{ig} + \varepsilon) - \ln \varepsilon}{\ln\left(\max_{g} x_{ig} + \varepsilon\right) - \ln \varepsilon}, \qquad \varepsilon = 10^{-6} \qquad (3)$$

Gen brez štetij ima po enačbi 3 vrednost 0, najbolj izražen gen v celici pa vrednost 1. Ker je $\lvert\ln \varepsilon\rvert$ velik v primerjavi z logaritmi štetij, so vrednosti vseh izraženih genov razmeroma blizu 1, zato normalizirana matrika v veliki meri zapiše, ali je gen v celici izražen. Normalizirano matriko smo z analizo glavnih komponent (angl. *principal component analysis*, PCA) iz paketa scikit-learn [@pedregosa2011] skrčili na 100 komponent, kot v izvirnem modelu [@he2026trim]. PCA smo izračunali na vseh celicah, ker ne uporablja oznak ekspanzije. Za prikaz celic v dveh razsežnostih smo uporabili metodo UMAP [@mcinnes2018]. Surova štetja smo shranili posebej, ker jih potrebuje scFEA (razdelek 3.4).

Klon smo definirali kot skupino celic z enakim zaporedjem CDR3β. Za vsak klon smo v vsakem kvadrantu prešteli celice. Enako kot v izvirni kodi modela TRIM velikosti klona nismo ločevali po bolnikih, zato celice redkih klonov z enakim zaporedjem CDR3β pri več bolnikih štejemo skupaj [@he2026trim]. Klon je ekspandiran, če ima v krvi po zdravljenju več celic kot v krvi pred zdravljenjem; vsaka celica dobi oznako svojega klona. Oznako smo prevzeli iz izvirne kode, ki primerja surova števila celic. V članku avtorji ekspanzijo opišejo s števili, normaliziranimi s številom vseh celic ob vsakem času [@he2026trim], zato smo v naknadni analizi preizkusili tudi takšno oznako (razdelek 3.9).

## Vložitev zaporedij TCR

Zaporedja CDR3β smo v vektorje pretvorili (vložili) s konvolucijskim avtokodirnikom iz izvirne kode modela TRIM [@he2026trim]. Zaporedja smo podaljšali na dolžino najdaljšega in vsako aminokislino zapisali s petimi Atchleyjevimi faktorji, ki povzemajo njene fizikalno-kemijske lastnosti [@atchley2005]. Kodirnik s tremi enodimenzionalnimi konvolucijskimi plastmi in polno povezano plastjo zaporedje preslika v vektor s 100 razsežnostmi, dekodirnik z obratno zgradbo pa iz vektorja rekonstruira Atchleyjeve faktorje. Avtokodirnik smo učili na enoličnih zaporedjih CDR3β, tako da je vsak klon prispeval en primer, s povprečno kvadratno napako rekonstrukcije in majhno kaznijo za velikost vložitev (tabela P2 v Prilogi). Kot v izvirni kodi smo vložitve pomnožili z 10, kadar je bila manj kot polovica njihovih razsežnosti po absolutni vrednosti večja od 1. Vsaka celica je dobila vložitev svojega klona.

## Ocenjevanje metabolnih fluksov s scFEA

Metabolne flukse smo ocenili z metodo scFEA (razdelek 1.7) [@alghamdi2021] in njenim modelom človeškega metabolizma s 168 metabolnimi moduli, povezanimi prek 70 presnovkov (datoteki `module_gene_m168.csv` in `cmMat_c70_m168.csv` iz repozitorija scFEA). Grafovska nevronska mreža se uči brez oznak, z minimizacijo vsote štirih členov [@alghamdi2021]: neravnovesja presnovkov med moduli, kazni za negativne flukse, odstopanja vsote fluksov celice od skupnega izražanja genov v celici in razlike med fluksi modulov celice ter povprečnim izražanjem genov teh modulov.

Na vhod smo dali surova štetja 658 od 663 genov modulov, ki so bili v naših podatkih. scFEA vhodno matriko logaritmira ($\log_2(x+1)$), kadar njena največja vrednost preseže 50. Imputacije manjkajočih vrednosti nismo uporabili, ostale nastavitve pa so bile privzete (tabela P2).

Izvorno kodo scFEA smo prilagodili novejšim različicam knjižnic in počasno zanko za sestavo vhoda mreže zamenjali z enakovrednim vektoriziranim izračunom; ujemanje smo preverili na manjšem vzorcu. Vseh celic kljub temu ni bilo mogoče obdelati naenkrat, zato smo jih v vrstnem redu iz predobdelave razdelili na 14 kosov po največ 5000 celic. Vsak kos je dobil svojo naučeno mrežo, skupno izražanje celic, izračunano samo iz genov modulov, pa je scFEA normaliziral s povprečjem znotraj kosa.

<!-- ZVEZEK: fluksi_kosi_sestava.csv (kateri bolniki in kvadranti so v katerem kosu), fluksi_kosi_eta2.csv. -->

Celoten izračun smo izvedli trikrat (semena 0, 1 in 2) in vsak zagon uporabili pri učenju modela TRIM z istim semenom. Ker scFEA nastavi lastno fiksno seme, smo pričakovali, da se bodo zagoni razlikovali le zaradi nedeterminističnih operacij na grafičnem procesorju. Rezultat scFEA je za vsako celico vektor 168 fluksov. scFEA flukse oceni samo iz izražanja genov, brez meritev hitrosti reakcij, zato so relativni in odvisni od modela; ne razlagamo jih kot absolutne hitrosti reakcij. Njihove vrednosti so bile majhne (povprečje okoli 0,006), zato je bila v razvojnih zagonih napaka njihove rekonstrukcije v funkciji izgube modela TRIM (razdelek 3.5) zanemarljiva in model je flukse prezrl. Vsak modul $m$ smo zato standardizirali na vseh celicah:

$$\tilde{f}_{im} = \frac{f_{im} - \bar{f}_m}{s_m + 10^{-8}} \qquad (4)$$

kjer sta $\bar{f}_m$ in $s_m$ povprečje in standardni odklon fluksa modula $m$ prek vseh celic. Standardizirane flukse (enačba 4) smo v model dali neposredno, brez PCA.

## Model TRIM in njegove različice

Model TRIM [@he2026trim] se nauči opisati vsako celico z latentnim vektorjem (stisnjenim opisom celice), iz katerega lahko nato generira celice v drugih pogojih, na primer po zdravljenju. Iz izvirne kode smo ga prenesli v zvezke Jupyter za okolje Google Colab. Zgradbo, funkcijo izgube in hiperparametre smo ohranili (tabela P2 v Prilogi), razlike od izvirnika opisujemo v razdelkih 3.6 in 3.7.

Za vsako modaliteto ima model svoj kodirnik in dekodirnik. Oba sta polno povezani nevronski mreži s tremi skritimi plastmi (2048, 1024 in 512 nevronov pri kodirniku in obratno pri dekodirniku), normalizacijo paketa in aktivacijo LeakyReLU. Kodirnik na vhodu dobi podatke modalitete in tri vektorje pogojev (tkivo, čas vzorčenja in bolnik, vsak s 128 razsežnostmi), na izhodu pa povprečje $\mu$ in logaritem variance $\log\sigma^2$ latentne porazdelitve s 1024 razsežnostmi. Model izhode kodirnikov vseh modalitet povpreči v skupni $\mu$ in $\log\sigma^2$, iz te porazdelitve vzorči latentni vektor $z$ in ga skupaj z vektorji pogojev preda dekodirnikom. Vektorje pogojev se model nauči sproti: pomožna mreža stalnega vzorca 5120 učnih celic preslika v vektorje, ki jih model povpreči po celicah iz istega tkiva, z istim časom vzorčenja oziroma istega bolnika iz krvi pred zdravljenjem [@he2026trim].

Funkcija izgube je povprečje petih uteženih členov (uteži so v tabeli P2). Člen $\mathcal{L}_{\mathrm{KL}}$ je Kullback-Leiblerjeva divergenca latentne porazdelitve od standardne normalne, povprečena po celicah in razsežnostih, in ustreza drugemu členu enačbe 1. $\mathcal{L}_{\mathrm{RNA}}$ je povprečna kvadratna napaka rekonstrukcije 100 komponent PCA, $\mathcal{L}_{\mathrm{vl}}$ pa povprečni kvadrat vektorjev pogojev. Receptor celic T model obravnava z dvema kontrastnima členoma. Dekodirnik TCR za vsako celico $i$ v paketu vrne vektor $\hat{t}_i$, celici $i$ in $j$ pa pripadata klonoma $c_i$ in $c_j$:

$$\mathcal{L}_{\mathrm{ist}} = \mathop{\mathrm{povp}}\limits_{i \neq j,\; c_i = c_j} \lVert \hat{t}_i - \hat{t}_j \rVert_1, \qquad \mathcal{L}_{\mathrm{razl}} = \mathop{\mathrm{povp}}\limits_{c_i \neq c_j} \max\left(0,\; \delta - \lVert \hat{t}_i - \hat{t}_j \rVert_1\right) \qquad (5)$$

Prvi člen enačbe 5 celice istega klona približa, drugi pa celice različnih klonov razmakne vsaj na razdaljo $\delta$ [@he2026trim]. Model se torej ne uči rekonstruirati vložitve TCR, temveč razdalje med kloni.

Primerjali smo tri različice modela, ki se razlikujejo samo po vhodnih modalitetah (tabela 1). Različica RNA + TCR je izvirni model [@he2026trim]. V različici fluksi + TCR smo kodirnik in dekodirnik RNA zamenjali s kodirnikom in dekodirnikom enake zgradbe za 168 standardiziranih fluksov, člen $\mathcal{L}_{\mathrm{RNA}}$ pa z napako rekonstrukcije fluksov $\mathcal{L}_{\mathrm{fluks}}$. Različica RNA + fluksi + TCR ima tri kodirnike in tri dekodirnike. Model izhode vseh treh kodirnikov povpreči v skupno latentno porazdelitev, funkcija izgube pa ima šesti člen $\lambda_{\mathrm{fluks}}\,\mathcal{L}_{\mathrm{fluks}}$. Pri vektorjih pogojev model povpreči pomožne mreže vseh modalitet. Vse ostale nastavitve so enake (tabela P2).

Model smo učili 20 000 korakov s paketi po 4096 celic in optimizatorjem Adam [@kingma2014adam] v knjižnici PyTorch [@paszke2019]. Model smo shranili po 10 000 in 20 000 korakih ter za napovedi uporabili končnega.

Table: Tabela 1: Različice modela TRIM in vhodi izhodiščnih klasifikatorjev.

| Različica | Vhodne modalitete | Izhodiščni klasifikatorji na |
|---|---|---|
| RNA + TCR | RNA (100 komponent PCA), TCR (100) | RNA, TCR, RNA in TCR |
| fluksi + TCR | fluksi (168), TCR (100) | fluksi, TCR, fluksi in TCR |
| RNA + fluksi + TCR | RNA (100), fluksi (168), TCR (100) | RNA, TCR, RNA in TCR |

## Napoved ekspanzije s psevdokloni

Model napove ekspanzijo tako, da za celice bolnika iz krvi pred zdravljenjem generira njihove »različice« po zdravljenju in primerja, kako velike klone tvorijo generirani receptorji pred zdravljenjem in po njem. Napoved smo izračunali po postopku iz izvirne kode [@he2026trim]. Za vsako celico $i$ izpuščenega bolnika iz krvi pred zdravljenjem smo s kodirniki izračunali $\mu_i$ in $\sigma_i$. Iz te porazdelitve smo $K$-krat vzorčili latentni vektor, kjer je $K$ celi del količnika $10\,000/n$ ($n$ je število celic, najmanj 1). Vsak vzorec smo dekodirali v vseh štirih kvadrantih, tako da smo zamenjali vektorja tkiva in časa, vektor bolnika pa ohranili. Generirane vektorje TCR vsakega kvadranta smo združili v psevdoklone: naključno izbrana še nerazvrščena celica in vse nerazvrščene celice, ki so ji bližje od $\tau/2$ v razdalji L1, tvorijo nov psevdoklon. Postopek smo ponavljali, dokler niso bile razvrščene vse celice. Za vsako celico smo nato v vsakem kvadrantu prešteli velikost njenega psevdoklona in jo povprečili po vzorcih. Iz povprečnih velikosti $\bar{p}_i$ v krvi pred zdravljenjem in po njem smo izračunali oceno ekspanzije:

$$y_i = \begin{cases} \bar{p}_i^{\,\mathrm{po}} / \bar{p}_i^{\,\mathrm{pred}}, & \bar{p}_i^{\,\mathrm{pred}} \geq 2 \\ 0, & \text{sicer} \end{cases} \qquad (6)$$

Ocena 0 v enačbi 6 pomeni, da model celico pred zdravljenjem uvrsti v klon z eno samo celico, kar kaže na neekspandiran klon.

Prag $\tau$ smo umerili tako, da je število psevdoklonov med celicami izpuščenega bolnika iz krvi pred zdravljenjem čim bližje številu resničnih klonov med njimi. Začeli smo z vrednostjo $\tau = \delta = 10$, vektorje TCR, ki jih je model za celice izpuščenega bolnika generiral iz naključnih latentnih vektorjev in pogojev teh celic, združili v psevdoklone znotraj istega tkiva in prešteli psevdoklone. Če jih je bilo premalo, smo prag pomnožili z 0,95, če preveč, pa z 1,05. Postopek smo končali, ko se je število psevdoklonov od števila klonov razlikovalo za največ 2,5 %, ali po 50 ponovitvah. Izvirna koda prag umeri na učnih celicah vseh bolnikov hkrati [@he2026trim], za kar bi matrika razdalj med vsemi celicami pri naših podatkih zasedla več kot 100 GB pomnilnika, zato smo ga umerili samo na celicah izpuščenega bolnika. Umerjanje uporablja le število klonov v krvi pred zdravljenjem, ne pa podatkov o ekspanziji, lahko pa nekoliko vpliva na AUC.

## Shema vrednotenja in izhodiščni klasifikatorji

**Shema »izpusti enega bolnika«.** Model smo vsakič učili brez enega bolnika. Učna množica je obsegala vse celice ostalih 27 bolnikov in celice izpuščenega bolnika iz krvi pred zdravljenjem. Te celice model potrebuje za vektor bolnika in so izhodišče napovedi, oznak ekspanzije pa ne vsebujejo. Celice izpuščenega bolnika iz krvi po zdravljenju in iz tumorja model med učenjem ni videl. Pri bolniku 23 je imela učna množica 66 829 celic, izpuščenih pa je bilo 3101 celic. Napoved smo vrednotili na vseh celicah izpuščenega bolnika iz krvi pred zdravljenjem.

**Izbira izpuščenih bolnikov.** AUC je zanesljiva le, če ima izpuščeni bolnik dovolj pozitivnih primerov. Kot merilo smo zato postavili vsaj 20 celic v krvi pred zdravljenjem, katerih klon se po zdravljenju namnoži. Poskuse smo zaradi računske zahtevnosti izvedli za dva bolnika, ki merilo izpolnjujeta. Bolnik 24 je imel celice v vseh štirih kvadrantih in med vsemi bolniki največ celic v tumorju po zdravljenju, zato smo ga izbrali za prvi poskus. Bolnik 23 je imel med vsemi bolniki največ celic iz ekspandiranih klonov v krvi pred zdravljenjem.

**Izhodiščni klasifikatorji.** Kot v izvirnem delu [@he2026trim] smo napoved modela primerjali z dvema klasifikatorjema iz paketa scikit-learn [@pedregosa2011]: metodo $k$ najbližjih sosedov (angl. *k-nearest neighbours*, kNN; $k = 20$) in večplastnim perceptronom s skritimi plastmi s 400, 200 in 100 nevroni (NN). Klasifikatorja smo učili na celicah iz krvi vseh ostalih bolnikov (pred zdravljenjem in po njem) z oznako, ali je klon celice ekspandiran. Testirali smo ju na celicah izpuščenega bolnika iz krvi pred zdravljenjem, katerih zaporedje CDR3β se v učni množici ne pojavi, in pri tem odstranili podvojene vrstice. Vhodi so bili RNA, TCR ali oboje skupaj, pri različici fluksi + TCR pa fluksi namesto RNA (tabela 1). <!-- PREVERI: koda, ki je dala izhodiščne klasifikatorje na fluksih (knn_flux, NN_flux v zapisih), ni v repozitoriju. Shranjeni zvezek 04_expansion_eval_flux ima starejšo različico. Ali so bili fluksi tu standardizirani? -->

**Semena in ponovitve.** Vsako kombinacijo bolnika in različice smo učili s tremi semeni (0, 1 in 2). Seme določi začetne uteži, vrstni red celic v paketih in vzorčenje latentnih vektorjev (knjižnici NumPy in PyTorch). Nekatere operacije na grafičnem procesorju niso deterministične, zato se lahko zagona z istim semenom nekoliko razlikujeta. Izvedli smo 17 od 18 načrtovanih zagonov: pri bolniku 24 manjka zagon različice RNA + TCR s semenom 2. Zagon različice fluksi + TCR pri bolniku 23 s semenom 0 smo izvedli dvakrat. V glavni analizi smo uporabili prvi zagon, drugega pa le za oceno ponovljivosti pri istem semenu. Vsak zagon smo zapisali v datoteko JSON z resničnimi oznakami, ocenami $y_i$ in rezultati izhodiščnih klasifikatorjev.

## Statistična analiza

Glavno merilo uspešnosti je bila AUC (razdelek 1.5) [@hanley1982], izračunana iz resnične oznake in ocene $y_i$. Za opis napovedi pri enem pragu smo uporabili prag z največjim Youdenovim indeksom (razlika med deležem pravilno pozitivnih in lažno pozitivnih) [@youden1950]. Prag smo izbrali na istih celicah, na katerih smo ga ovrednotili.

Za vsakega bolnika in različico smo poročali AUC posameznih semen, njihovo povprečje in razpon ter AUC ansambla semen. Ansambel smo izračunali tako, da smo za vsako celico povprečili range njenih ocen $y_i$ v zagonih z različnimi semeni (vezi so dobile povprečni rang). Interval zaupanja za AUC smo izračunali na dva načina. DeLongova metoda [@delong1988] predpostavi neodvisne celice. Celice istega klona pa imajo isto oznako in podobne ocene, zato smo izračunali tudi interval s ponovnim vzorčenjem s ponavljanjem (angl. *bootstrap*) po klonih [@field2007]: 2000-krat smo s ponavljanjem izbrali toliko klonov, kot jih je imel bolnik, vzeli vse njihove celice in iz porazdelitve AUC odčitali 2,5. in 97,5. percentil. Pred delom smo predvideli DeLongov interval; ker se je izkazalo, da ga odvisnost celic istega klona močno zoži, smo za sklepe uporabili interval bootstrap po klonih, DeLongovega pa navajamo za primerjavo.

AUC modela TRIM je bila višja od AUC izhodiščnega klasifikatorja, če je bila spodnja meja njenega 95-% intervala zaupanja nad AUC najboljšega izhodiščnega klasifikatorja iz istega zagona. Izhodiščna klasifikatorja smo testirali le na celicah z zaporedji CDR3β, ki jih v učni množici ni bilo, zato smo AUC modela izračunali tudi na teh celicah. Napovedi izhodiščnih klasifikatorjev nismo shranili, zato parni test zanje ni bil mogoč. Različici smo primerjali na istih celicah pri istem bolniku in semenu z razliko v AUC in njenim 95-% intervalom ter vrednostjo *p* iz parnega bootstrapa po klonih (ista izbira klonov za obe oceni); za primerjavo smo izračunali tudi DeLongov test za korelirane krivulje ROC [@delong1988]. Za hipotezo H3 pred delom nismo določili natančnejšega pravila. Izboljšanje smo razumeli kot višjo AUC različice s fluksi, ki bi bila statistično značilna in istega predznaka pri vseh semenih in obeh bolnikih. Z enakima testoma smo primerjali tudi zagone iste različice z različnimi semeni. Teh 16 primerjav je opisnih, zato jih nismo popravili za večkratno testiranje. Podobnost razvrstitve celic med zagoni smo ocenili s Spearmanovim koeficientom korelacije ocen $y_i$ [@spearman1904].

Ujemanje fluksov med zagoni scFEA smo ocenili s Pearsonovim koeficientom korelacije [@pearson1895]. Delež variance fluksa modula, ki ga pojasni pripadnost kosu, bolniku ali kvadrantu, smo izračunali kot razmerje med vsoto kvadratov med skupinami in skupno vsoto kvadratov ($\eta^2$). Flukse posameznih modulov smo med ekspandiranimi in neekspandiranimi kloni v krvi pred zdravljenjem primerjali na ravni klonov: fluks klona je bil povprečje standardiziranih fluksov njegovih celic, kot enoto pa smo šteli par bolnika in zaporedja CDR3β. Klone smo primerjali s testom Manna in Whitneyja [@mann1947] pri vseh bolnikih, samo med kloni z eno celico in ločeno za celice T CD4+ in CD8+, vrednosti $p$ pa popravili za večkratno testiranje po Benjaminiju in Hochbergu [@benjamini1995] čez 168 modulov. Iz statistike $U$ smo izračunali tudi AUC posameznega modula, $U/(n_1 n_0)$, kjer sta $n_1$ in $n_0$ števili klonov v obeh skupinah.

## Velikost klona kot izhodišče (naknadna analiza)

Ker oznaka ekspanzije primerja število celic klona v krvi po zdravljenju in pred njim, smo po zaključenih zagonih preverili, koliko AUC modela TRIM pojasni že velikost klona pred zdravljenjem; te analize pred delom nismo načrtovali. Kot oceno za vsako celico smo vzeli število celic njenega klona v krvi pred zdravljenjem, prešteto prek vseh bolnikov kot pri oznaki. Za primerjavo smo izračunali tudi oceni iz števila celic klona samo pri izpuščenem bolniku in iz števila celic klona v tumorju pred zdravljenjem. Ansambel vseh zagonov je združil vseh 9 zagonov pri bolniku 23 oziroma 8 pri bolniku 24. AUC smo izračunali na vseh celicah, ločeno na celicah iz klonov z eno celico in z vsaj dvema celicama v krvi pred zdravljenjem ter na ravni klonov, kjer je vsak klon štel enkrat s povprečno oceno svojih celic. Na ravni celic smo uporabili intervale in parne primerjave bootstrap po klonih (razdelek 3.8), na ravni klonov pa DeLongov interval in parni DeLongov test [@delong1988]. Da bi preverili, ali model razlikuje med kloni enake velikosti, smo celice razvrstili po velikosti klona in celice z enako velikostjo razvrstili še po oceni modela.

Preizkusili smo tudi oznako iz deležev: klon je bil ekspandiran, če je bil njegov delež med celicami bolnika v krvi po zdravljenju večji kot pred zdravljenjem. Povezavo med deležem ekspandiranih celic in razmerjem med številom celic v krvi po zdravljenju in pred njim smo pri 24 bolnikih z obema vzorcema krvi ocenili s Spearmanovim koeficientom korelacije [@spearman1904].

## Programska in strojna oprema

Analizo smo izvedli v programskem jeziku Python s knjižnicami scanpy [@wolf2018], PyTorch [@paszke2019], scikit-learn [@pedregosa2011], umap-learn [@mcinnes2018], NumPy, pandas in SciPy ter z izvorno kodo scFEA [@alghamdi2021]. Modele smo učili v okolju Google Colab z grafičnim procesorjem NVIDIA A100 z 80 GB pomnilnika in razširjenim delovnim pomnilnikom (okoli 83 GB). Verzije knjižnic vpišemo po zagonu statističnega zvezka. <!-- ZVEZEK: verzije.json. Opomba: verzije so iz okolja Colab ob zagonu zvezka; poskuse smo izvedli julija in avgusta 2026. Dodaj čas učenja enega zagona na A100 (lab_diary.md), če je zapisan. -->

Vsa koda, zapisi zagonov in laboratorijski dnevnik so v repozitoriju <https://github.com/azi-kos/TRIM-Flux> (veja `flux-extension`). Statistične analize in slike smo izdelali s skriptami in zvezki v mapi `diploma/skripte`. <!-- Ob oddaji dodaj oznako zadnje potrditve (commit). -->
