<!-- Materiali in metode (Navodila 2026 + Biokemija v1.2): dovolj natančno za ponovitev, logični vrstni red,
     vsaka metoda iz literature z virom, verzije programske opreme in strojna oprema.
     Vse številke so iz kode (notebooks*/), results/*/args.txt, results/*/heldout_meta.json, results/logs/logs/*.json
     in lab_diary.md. Številke z oznako »ZVEZEK« pridejo iz diploma/skripte/colab_statistika.ipynb. -->

# Materiali in metode

Potek dela od surovih podatkov do vrednotenja prikazuje slika 4. Vso analizo smo izvedli v zvezkih Jupyter v okolju Google Colab; koda je v repozitoriju TRIM-Flux (razdelek 3.8).

<!-- SLIKA 4 (lastna, shema): GEO → kontrola kakovosti in normalizacija → (RNA: PCA 100) / (TCR: CDR3β → Atchley → CNN, 100) /
     (surova štetja → scFEA, 168 modulov → z) → TRIM (tri različice) → psevdokloni → AUC na izpuščenem bolniku. -->

**[Slika 4: potek dela od surovih podatkov do vrednotenja napovedi. Narišem jo po pregledu besedila.]**

## Podatki

Uporabili smo javno dostopne podatke raziskave Luome in sodelavcev [@luoma2022] iz zbirke Gene Expression Omnibus (GEO) z oznako GSE200996 [@barrett2013]. Raziskava je zajela bolnike s HNSCC, ki so pred operacijo prejeli neoadjuvantno imunoterapijo [@luoma2022]. Iz zbirke smo prenesli matrike štetij prepisov (datoteke `.h5`), anotacije kontigov TCR (datoteke `filtered_contig_annotations`) in metapodatke celic T CD4+ in CD8+ iz krvi in tumorja. Metapodatki za vsako celico navajajo bolnika, tkivo, podtip in oznako vzorca.

Vzorce z oznakama B1 (kri) in Pre-Tx (tumor) smo obravnavali kot vzorce pred zdravljenjem, vzorce z oznakama B2 (kri) in Post-Tx (tumor) pa kot vzorce po zdravljenju. Vzorcev z oznako B3 nismo uporabili, ker zanje zbirka ne vsebuje matrik štetij. Kombinacijo tkiva (kri ali tumor) in časa vzorčenja (pred zdravljenjem ali po njem) imenujemo kvadrant. Bolnike smo oštevilčili od 0 do 27 po abecednem vrstnem redu njihovih oznak v zbirki GEO. V nadaljevanju jih navajamo s temi indeksi, oznake GEO pa so v tabeli 1.

<!-- ZVEZEK: tabela1_celice_po_bolnikih.csv (bolnik_indeks, bolnik_id, kvadranti, CD4/CD8, kloni, ekspandirane celice).
     Število bolnikov (28) preveri v povzetek_podatkov.json. -->

**[Tabela 1: število celic po bolnikih in kvadrantih po kontroli kakovosti. Vstavim jo po zagonu statističnega zvezka.]**

## Predobdelava podatkov

Podatke smo predobdelali s paketom scanpy [@wolf2018] po splošnih priporočilih za analizo podatkov scRNA-seq [@luecken2019]. Datoteke `.h5` vsebujejo surove matrike štetij, v katerih so poleg celic tudi prazne mikrokapljice in poškodovane celice. Kontrolo kakovosti smo zato izvedli za vsako datoteko posebej v treh korakih:

1. Obdržali smo celice z vsaj 200 izraženimi geni. Prag smo izbrali na podlagi porazdelitve števila izraženih genov na celico (razdelek 4.1).
2. Obdržali smo samo celice, ki so v metapodatkih avtorjev [@luoma2022]. Avtorji so vzorce več bolnikov sekvencirali skupaj, zato se ista črtna koda lahko pojavi pri več bolnikih. Celico smo zato vedno označili s parom (črtna koda, bolnik) in ne s črtno kodo samo.
3. Odstranili smo celice, pri katerih je bil delež štetij mitohondrijskih genov (imena z začetkom MT-) 20 % ali več, ker visok delež kaže na poškodovano celico [@luecken2019].

Iz anotacij kontigov TCR smo obdržali verige β (TRB), ki so bile produktivne, visoko zanesljive in so imele določeno zaporedje CDR3. Če je imela celica več takšnih kontigov, smo obdržali prvega. Izločili smo celice brez zaporedja CDR3β. Po vseh korakih je ostalo 69 930 celic s 36 601 genom (tabela 1).

<!-- ZVEZEK: qc_po_datotekah.csv za število celic, ki jih odstrani vsak korak. -->

Par (črtna koda, bolnik) ne loči tkiva in časa vzorčenja istega bolnika. V metapodatkih kaže 102 para na kri in tumor hkrati, 141 parov pa na več časov vzorčenja. V takih primerih smo obdržali prvo oznako, zato ima lahko napačno oznako tkiva ali časa največ okoli 0,3 % celic.

Štetja smo normalizirali enako kot avtorji modela TRIM [@he2026trim]. Za vsako celico $i$ smo štetja $x_{ig}$ logaritmirali in jih preslikali na interval od 0 do 1:

$$\tilde{x}_{ig} = \frac{\ln(x_{ig} + \varepsilon) - \ln \varepsilon}{\ln\left(\max_{g} x_{ig} + \varepsilon\right) - \ln \varepsilon}, \qquad \varepsilon = 10^{-6} \qquad (3)$$

Gen brez štetij ima po enačbi 3 vrednost 0, najbolj izražen gen v celici pa vrednost 1. Ker je $\lvert\ln \varepsilon\rvert$ velik v primerjavi z logaritmi štetij, so vrednosti vseh izraženih genov razmeroma blizu 1, zato normalizirana matrika v veliki meri zapiše, ali je gen v celici izražen. Normalizirano matriko smo z analizo glavnih komponent (angl. *principal component analysis*, PCA) iz paketa scikit-learn [@pedregosa2011] skrčili na 100 komponent, kot v izvirnem modelu [@he2026trim]. PCA smo izračunali na vseh celicah, ker ne uporablja oznak ekspanzije. Za prikaz celic v dveh razsežnostih smo uporabili metodo UMAP [@mcinnes2018]. Surova štetja smo shranili posebej, ker jih potrebuje scFEA (razdelek 3.5).

Klon smo definirali kot skupino celic z enakim zaporedjem CDR3β. Za vsak klon smo v vsakem kvadrantu prešteli celice. Enako kot v izvirni kodi modela TRIM velikosti klona nismo ločevali po bolnikih, zato celice redkih klonov z enakim zaporedjem CDR3β pri več bolnikih štejemo skupaj [@he2026trim].

## Vložitev zaporedij TCR

Zaporedja CDR3β smo v vektorje pretvorili s konvolucijskim avtokodirnikom iz izvirne kode modela TRIM [@he2026trim]. Zaporedja smo s presledki podaljšali na dolžino najdaljšega zaporedja in vsako aminokislino zapisali s petimi Atchleyjevimi faktorji, ki povzemajo fizikalno-kemijske lastnosti aminokislin [@atchley2005]. Kodirnik ima tri enodimenzionalne konvolucijske plasti (velikost jedra 3, korak 2, 512, 256 in 128 kanalov, normalizacija paketa in aktivacija LeakyReLU) in polno povezano plast, ki zaporedje preslika v vektor s 100 razsežnostmi. Dekodirnik z obratno zgradbo iz vektorja rekonstruira Atchleyjeve faktorje.

Avtokodirnik smo učili na enoličnih zaporedjih CDR3β, tako da je vsak klon prispeval en primer. Funkcija izgube je bila povprečna kvadratna napaka rekonstrukcije, ki smo ji prišteli 0,001-kratnik povprečnega kvadrata vložitve, da vložitve ne bi neomejeno rasle. Učenje je trajalo 50 100 korakov s paketi po 1024 zaporedij in optimizatorjem Adam [@kingma2014adam] s hitrostjo učenja 0,001. Kot v izvirni kodi smo vložitve pomnožili z 10, kadar je bila manj kot polovica njihovih razsežnosti po absolutni vrednosti večja od 1. Vsaka celica je dobila vložitev svojega klona.

## Model TRIM

Model TRIM [@he2026trim] smo iz izvirne kode prenesli v zvezke Jupyter za okolje Google Colab. Zgradbo, funkcijo izgube in hiperparametre smo ohranili (tabela 2), razlike od izvirnika opisujemo v razdelku 3.7.

Za vsako modaliteto ima model svoj kodirnik in dekodirnik. Oba sta polno povezani nevronski mreži s tremi skritimi plastmi (2048, 1024 in 512 nevronov pri kodirniku in obratno pri dekodirniku), normalizacijo paketa in aktivacijo LeakyReLU. Kodirnik na vhodu dobi podatke modalitete in tri vektorje pogojev (tkivo, čas vzorčenja in bolnik, vsak s 128 razsežnostmi), na izhodu pa povprečje $\mu$ in logaritem variance $\log\sigma^2$ latentne porazdelitve s 1024 razsežnostmi. Model izhode kodirnikov vseh modalitet povpreči v skupni $\mu$ in $\log\sigma^2$, iz te porazdelitve vzorči latentni vektor $z$ in ga skupaj z vektorji pogojev preda dekodirnikom.

Vektorje pogojev se model nauči sproti. Na začetku učenja smo naključno izbrali 5120 učnih celic in jih obdržali ves čas učenja. V vsakem koraku pomožna mreža vsake modalitete te celice preslika v vektorje s 128 razsežnostmi, model pa vektorje modalitet povpreči. Vektor pogoja za tkivo je povprečje teh vektorjev po celicah iz istega tkiva, vektor za čas je povprečje po celicah z istim časom vzorčenja, vektor bolnika pa povprečje po celicah tega bolnika iz krvi pred zdravljenjem.

Funkcija izgube je povprečje petih uteženih členov:

$$\mathcal{L} = \frac{1}{5}\left(\lambda_{\mathrm{KL}}\,\mathcal{L}_{\mathrm{KL}} + \lambda_{\mathrm{RNA}}\,\mathcal{L}_{\mathrm{RNA}} + \lambda_{\mathrm{TCR}}\,\mathcal{L}_{\mathrm{ist}} + \lambda_{\mathrm{TCR}}\,\mathcal{L}_{\mathrm{razl}} + \lambda_{\mathrm{vl}}\,\mathcal{L}_{\mathrm{vl}}\right) \qquad (4)$$

Člen $\mathcal{L}_{\mathrm{KL}} = -\frac{1}{nd}\sum_{i}\sum_{k}\left(1 + \log\sigma_{ik}^2 - \sigma_{ik}^2 - \mu_{ik}^2\right)$ približuje latentno porazdelitev standardni normalni ($n$ je število celic v paketu, $d$ razsežnost latentnega prostora), kar ustreza drugemu členu enačbe 1. $\mathcal{L}_{\mathrm{RNA}}$ je povprečna kvadratna napaka rekonstrukcije 100 komponent PCA, $\mathcal{L}_{\mathrm{vl}}$ pa povprečni kvadrat vektorjev pogojev. Receptor celic T model obravnava s kontrastno izgubo. Dekodirnik TCR za vsako celico $i$ v paketu vrne vektor $\hat{t}_i$, celici $i$ in $j$ pa pripadata klonoma $c_i$ in $c_j$:

$$\mathcal{L}_{\mathrm{ist}} = \mathop{\mathrm{povp}}\limits_{i \neq j,\; c_i = c_j} \lVert \hat{t}_i - \hat{t}_j \rVert_1, \qquad \mathcal{L}_{\mathrm{razl}} = \mathop{\mathrm{povp}}\limits_{c_i \neq c_j} \max\left(0,\; \delta - \lVert \hat{t}_i - \hat{t}_j \rVert_1\right) \qquad (5)$$

Prvi člen enačbe 5 celice istega klona približa, drugi pa celice različnih klonov razmakne vsaj na razdaljo $\delta$ [@he2026trim]. Model se torej ne uči rekonstruirati vložitve TCR, temveč razdalje med kloni.

Model smo učili 20 000 korakov s paketi po 4096 celic in optimizatorjem Adam [@kingma2014adam] (tabela 2) v knjižnici PyTorch [@paszke2019]. Model smo shranili po 10 000 in 20 000 korakih ter za napovedi uporabili končnega.

Table: Tabela 2: Hiperparametri modela TRIM, enaki za vse različice in bolnike (iz datotek `args.txt`).

| Hiperparameter | Vrednost |
|---|---:|
| število korakov učenja | 20 000 |
| velikost paketa | 4096 |
| hitrost učenja (Adam) | 0,001 |
| razsežnost latentnega prostora $d$ | 1024 |
| razsežnost vektorja pogoja | 128 |
| število celic za vektorje pogojev | 5120 |
| utež $\lambda_{\mathrm{KL}}$ | 15 |
| uteži $\lambda_{\mathrm{RNA}}$, $\lambda_{\mathrm{fluks}}$, $\lambda_{\mathrm{TCR}}$ | 1 |
| utež $\lambda_{\mathrm{vl}}$ | 10 |
| rob kontrastne izgube $\delta$ | 10 |

## Ocenjevanje metabolnih fluksov s scFEA

Metabolne flukse smo ocenili z metodo scFEA [@alghamdi2021]. Uporabili smo njen model človeškega metabolizma s 168 metabolnimi moduli, ki so med seboj povezani prek 70 presnovkov (datoteki `module_gene_m168.csv` in `cmMat_c70_m168.csv` iz repozitorija scFEA). Za vsak modul je znan seznam genov, katerih produkti katalizirajo njegove reakcije. scFEA za vsak nabor podatkov nauči grafovsko nevronsko mrežo, ki iz izražanja genov modula oceni fluks skozi modul v posamezni celici. Mreža se uči brez oznak, z minimizacijo vsote štirih členov [@alghamdi2021]: neravnotežja presnovkov med moduli, kazni za negativne flukse, odstopanja vsote fluksov celice od skupnega izražanja genov v celici in razlike med fluksi modulov celice ter povprečnim izražanjem genov teh modulov.

Na vhod smo dali surova štetja genov, ki so v seznamih genov modulov <!-- ZVEZEK: scfea_geni.json, nasih_genov_v_modulih -->. scFEA vhodno matriko logaritmira ($\log_2(x+1)$), kadar njena največja vrednost preseže 50. Imputacije manjkajočih vrednosti nismo uporabili. Ostale nastavitve so bile privzete: 100 epoh, hitrost učenja 0,008 in optimizator Adam, pri čemer scFEA v vsaki epohi uporabi vse celice naenkrat. scFEA v izvorni kodi nastavi fiksno seme generatorja naključnih števil (16).

Izvorna koda scFEA od leta 2021 ni posodobljena, zato smo jo za novejše različice knjižnic popravili. Zanko, ki iz izražanja genov sestavi vhod mreže in je pri skoraj 70 000 celicah trajala več ur in porabila preveč pomnilnika, smo zamenjali z enakovrednim vektoriziranim izračunom in ujemanje rezultatov preverili na manjšem vzorcu. Tenzorje z grafičnega procesorja smo pred zapisom prenesli v glavni pomnilnik, zastarele funkcije knjižnice pandas pa zamenjali z novejšimi. Kljub temu vseh celic ni bilo mogoče obdelati naenkrat, zato smo jih v vrstnem redu iz predobdelave razdelili na 14 kosov po največ 5000 celic in scFEA pognali za vsak kos posebej. Vsak kos je tako dobil svojo naučeno mrežo, skupno izražanje celic pa je scFEA normaliziral s povprečjem znotraj kosa. Ker smo na vhod dali samo gene iz modulov, je scFEA skupno izražanje celice izračunal iz teh genov in ne iz vseh genov. Posledice obeh odločitev opisujemo v Razpravi.

<!-- ZVEZEK: fluksi_kosi_sestava.csv (kateri bolniki in kvadranti so v katerem kosu), fluksi_kosi_eta2.csv. -->

Celoten izračun smo izvedli trikrat in zagone označili s semeni 0, 1 in 2. Vsak zagon smo uporabili pri učenju modela TRIM z istim semenom. Ker scFEA nastavi lastno fiksno seme, smo pričakovali, da se bodo zagoni razlikovali le zaradi nedeterminističnih operacij na grafičnem procesorju. Ujemanje med zagoni smo preverili s Pearsonovim korelacijskim koeficientom (razdelek 4.4).

Rezultat scFEA je za vsako celico vektor 168 fluksov. scFEA flukse oceni samo iz izražanja genov, brez meritev hitrosti reakcij, zato so relativni in odvisni od modela; ne razlagamo jih kot absolutne hitrosti reakcij. Njihove vrednosti so bile majhne (povprečje okoli 0,006), zato je bila napaka rekonstrukcije fluksov v funkciji izgube zanemarljiva v primerjavi s členi za TCR in model flukse prezrl. Vsak modul $m$ smo zato standardizirali na vseh celicah:

$$\tilde{f}_{im} = \frac{f_{im} - \bar{f}_m}{s_m + 10^{-8}} \qquad (6)$$

kjer sta $\bar{f}_m$ in $s_m$ povprečje in standardni odklon fluksa modula $m$ prek vseh celic. Standardizirane flukse smo v model dali neposredno, brez PCA.

## Različice modela

Primerjali smo tri različice modela TRIM, ki se razlikujejo samo po vhodnih modalitetah (tabela 3). Različica RNA + TCR je izvirni model [@he2026trim]. V različici fluksi + TCR smo kodirnik in dekodirnik RNA zamenjali s kodirnikom in dekodirnikom enake zgradbe za 168 standardiziranih fluksov, napako rekonstrukcije RNA v enačbi 4 pa z napako rekonstrukcije fluksov $\mathcal{L}_{\mathrm{fluks}}$. Različica RNA + fluksi + TCR ima tri kodirnike in tri dekodirnike. Model izhode vseh treh kodirnikov povpreči v skupno latentno porazdelitev, funkcija izgube pa je povprečje šestih členov: enačbi 4 doda člen $\lambda_{\mathrm{fluks}}\,\mathcal{L}_{\mathrm{fluks}}$. Pri vektorjih pogojev model povpreči pomožne mreže vseh modalitet. Vse ostale nastavitve so enake (tabela 2).

Table: Tabela 3: Različice modela TRIM in vhodi izhodiščnih klasifikatorjev.

| Različica | Oznaka v kodi | Vhodne modalitete | Izhodiščni klasifikatorji na |
|---|---|---|---|
| RNA + TCR | `rna_tcr` | RNA (100 komponent PCA), TCR (100) | RNA, TCR, RNA in TCR |
| fluksi + TCR | `flux_tcr` | fluksi (168), TCR (100) | fluksi, TCR, fluksi in TCR |
| RNA + fluksi + TCR | `rna_flux_tcr` | RNA (100), fluksi (168), TCR (100) | RNA, TCR, RNA in TCR |

## Vrednotenje napovedi klonske ekspanzije

**Shema »izpusti enega bolnika«.** Model smo vsakič učili brez enega bolnika. Učna množica je obsegala vse celice ostalih 27 bolnikov in celice izpuščenega bolnika iz krvi pred zdravljenjem. Te celice model potrebuje za vektor bolnika in so izhodišče napovedi, oznak ekspanzije pa ne vsebujejo. Celice izpuščenega bolnika iz krvi po zdravljenju in iz tumorja model med učenjem ni videl. Pri bolniku 23 je imela učna množica 66 829 celic, izpuščenih pa je bilo 3101 celic.

**Izbira izpuščenih bolnikov.** AUC je zanesljiva le, če ima izpuščeni bolnik dovolj pozitivnih primerov. Kot merilo smo zato postavili vsaj 20 celic v krvi pred zdravljenjem, katerih klon se po zdravljenju razširi. Poskuse smo zaradi računske zahtevnosti izvedli za dva bolnika, ki merilo izpolnjujeta. Bolnik 24 je imel celice v vseh štirih kvadrantih in med vsemi bolniki največ celic v tumorju po zdravljenju, zato smo ga izbrali za prvi poskus. Bolnik 23 je bil med bolniki z največ ekspandiranimi celicami v krvi pred zdravljenjem. <!-- ZVEZEK: preveri v tabeli 1, ali ima bolnik 23 največ ekspandiranih celic. -->

**Resnična ekspanzija.** Klon je ekspandiran, če ima v krvi po zdravljenju več celic kot v krvi pred zdravljenjem. Napoved smo vrednotili na vseh celicah izpuščenega bolnika iz krvi pred zdravljenjem; vsaka celica je dobila oznako svojega klona. Pri bolniku 23 je bilo takih celic 1159, od tega 364 (31,4 %) iz ekspandiranih klonov, pri bolniku 24 pa 523, od tega 23 (4,4 %).

**Napoved modela TRIM.** Napoved smo izračunali po postopku iz izvirne kode [@he2026trim]. Za vsako celico $i$ izpuščenega bolnika iz krvi pred zdravljenjem smo s kodirniki izračunali $\mu_i$ in $\sigma_i$. Iz te porazdelitve smo $S$-krat vzorčili latentni vektor, kjer je $S$ celi del količnika $10\,000/n$ ($n$ je število celic, najmanj 1). Vsak vzorec smo dekodirali v vseh štirih kvadrantih, tako da smo zamenjali vektorja tkiva in časa, vektor bolnika pa ohranili. Generirane vektorje TCR vsakega kvadranta smo združili v psevdoklone: naključno izbrana še nerazvrščena celica in vse nerazvrščene celice, ki so ji bližje od $\tau/2$ v razdalji L1, tvorijo nov psevdoklon. Postopek smo ponavljali, dokler niso bile razvrščene vse celice. Za vsako celico smo nato v vsakem kvadrantu prešteli velikost njenega psevdoklona in jo povprečili po vzorcih. Iz povprečnih velikosti $\bar{p}_i$ v krvi pred zdravljenjem in po njem smo izračunali oceno ekspanzije:

$$y_i = \begin{cases} \bar{p}_i^{\,\mathrm{po}} / \bar{p}_i^{\,\mathrm{pred}}, & \bar{p}_i^{\,\mathrm{pred}} \geq 2 \\ 0, & \text{sicer} \end{cases} \qquad (7)$$

Ocena 0 v enačbi 7 pomeni, da model celico pred zdravljenjem uvrsti v klon z eno samo celico, kar kaže na neekspandiran klon.

**Umerjanje praga psevdoklonov.** Prag $\tau$ smo umerili tako, da je število psevdoklonov med celicami izpuščenega bolnika iz krvi pred zdravljenjem čim bližje številu resničnih klonov med njimi. Začeli smo z vrednostjo $\tau = \delta = 10$, generirane vektorje TCR celic izpuščenega bolnika združili v psevdoklone znotraj istega tkiva in prešteli psevdoklone. Če jih je bilo premalo, smo prag pomnožili z 0,95, če preveč, pa z 1,05. Postopek smo končali, ko se je število psevdoklonov od števila klonov razlikovalo za največ 2,5 %, ali po 50 ponovitvah. Izvirna koda prag umeri na učnih celicah vseh bolnikov hkrati [@he2026trim]. Za to potrebuje matriko razdalj med vsemi celicami, ki bi pri naših podatkih zasedla več kot 100 GB pomnilnika. Prag smo zato umerili samo na celicah izpuščenega bolnika. Umerjanje uporablja le število klonov v krvi pred zdravljenjem, ne pa podatkov o ekspanziji, lahko pa nekoliko vpliva na AUC.

**Izhodiščni klasifikatorji.** Kot v izvirnem delu [@he2026trim] smo napoved modela primerjali z dvema klasifikatorjema iz paketa scikit-learn [@pedregosa2011]: metodo $k$ najbližjih sosedov (angl. *k-nearest neighbours*, kNN; $k = 20$) in večplastnim perceptronom s skritimi plastmi s 400, 200 in 100 nevroni (NN). Klasifikatorja smo učili na celicah iz krvi vseh ostalih bolnikov (pred zdravljenjem in po njem) z oznako, ali je klon celice ekspandiran. Testirali smo ju na celicah izpuščenega bolnika iz krvi pred zdravljenjem, katerih zaporedje CDR3β se v učni množici ne pojavi, in pri tem odstranili podvojene vrstice. Vhodi so bili RNA, TCR ali oboje skupaj, pri različici fluksi + TCR pa fluksi namesto RNA (tabela 3).

**Merila uspešnosti.** Glavno merilo je bila ploščina pod krivuljo ROC, izračunana iz resnične oznake in ocene $y_i$. AUC je verjetnost, da model naključno izbrani celici iz ekspandiranega klona pripiše višjo oceno kot naključno izbrani celici iz neekspandiranega klona [@hanley1982]; vrednost 0,5 ustreza naključnemu ugibanju. Za matrike zamenjav smo oceno $y_i$ pretvorili v napoved pri treh pragovih: $y_i > 0$, pri pragu z največjim Youdenovim indeksom (razlika med deležem pravilno pozitivnih in lažno pozitivnih) [@youden1950] in pri pragu, ki kot ekspandirane označi enak delež celic, kot je resnično ekspandiranih.

**Semena in ponovitve.** Vsako kombinacijo bolnika in različice smo učili s tremi semeni (0, 1 in 2). Seme določi začetne uteži, vrstni red celic v paketih in vzorčenje latentnih vektorjev (knjižnici NumPy in PyTorch). Nekatere operacije na grafičnem procesorju niso deterministične, zato se lahko zagona z istim semenom nekoliko razlikujeta. Izvedli smo 17 od 18 načrtovanih zagonov: pri bolniku 24 manjka zagon različice RNA + TCR s semenom 2. Zagon različice fluksi + TCR pri bolniku 23 s semenom 0 smo izvedli dvakrat. V glavni analizi smo uporabili prvi zagon, drugega pa le za oceno ponovljivosti pri istem semenu. Vsak zagon smo zapisali v datoteko JSON z resničnimi oznakami, ocenami $y_i$ in rezultati izhodiščnih klasifikatorjev.

**Statistična analiza.** Za vsakega bolnika in različico smo poročali AUC posameznih semen ter njihovo povprečje in razpon. Interval zaupanja za AUC smo izračunali po DeLongovi metodi [@delong1988]. Različici smo primerjali na istih celicah pri istem bolniku in semenu z DeLongovim testom za korelirane krivulje ROC [@delong1988]. Z istim testom smo primerjali tudi zagone iste različice z različnimi semeni, podobnost razvrstitve celic med zagoni pa smo ocenili s Spearmanovim koeficientom korelacije ocen $y_i$ [@spearman1904]. Izhodiščne klasifikatorje smo testirali na drugi množici celic kot model TRIM, zato parni test zanje ni bil mogoč. AUC modela TRIM je bila višja od AUC izhodiščnega klasifikatorja, če je bila spodnja meja njenega 95-% intervala zaupanja nad AUC najboljšega izhodiščnega klasifikatorja iz istega zagona. Razlike v fluksih posameznih modulov med celicami iz ekspandiranih in neekspandiranih klonov v krvi pred zdravljenjem smo ocenili s testom Manna in Whitneyja [@mann1947] in vrednosti $p$ popravili za večkratno testiranje po Benjaminiju in Hochbergu [@benjamini1995]. Iz statistike $U$ smo izračunali tudi AUC posameznega modula, $U/(n_1 n_0)$, kjer sta $n_1$ in $n_0$ števili celic v obeh skupinah. Številke iz zapisov zagonov smo izračunali s skripto `rezultati_napovedi.py`, analize fluksov pa s statističnim zvezkom; oba sta v mapi `diploma/skripte` repozitorija.

## Programska in strojna oprema

Analizo smo izvedli v programskem jeziku Python s knjižnicami scanpy [@wolf2018], PyTorch [@paszke2019], scikit-learn [@pedregosa2011], umap-learn [@mcinnes2018], NumPy, pandas in SciPy ter z izvorno kodo scFEA [@alghamdi2021]. Modele smo učili v okolju Google Colab z grafičnim procesorjem NVIDIA A100 z 80 GB pomnilnika in razširjenim delovnim pomnilnikom (okoli 83 GB). Verzije knjižnic vpišemo po zagonu statističnega zvezka. <!-- ZVEZEK: verzije.json. Opomba: verzije so iz okolja Colab ob zagonu zvezka; poskuse smo izvedli julija in avgusta 2026. -->

Vsa koda, zapisi zagonov in laboratorijski dnevnik so v repozitoriju <https://github.com/azi-kos/TRIM-Flux> (veja `flux-extension`). <!-- Ob oddaji dodaj oznako zadnje potrditve (commit). -->
