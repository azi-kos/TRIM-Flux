<!-- ~1–2 strani. Pregled, ocena hipotez (potrjena/zavrnjena), težave, nadaljnje delo. -->

# Zaključek

Hipotezo H1 smo potrdili: model TRIM z RNA in TCR je ekspanzijo klonov v krvi izpuščenega bolnika v štirih od petih zagonov napovedal bolje kot izhodiščni klasifikatorji, tudi na istih celicah. Hipoteze H2 nismo potrdili, saj je model s fluksi in TCR izhodiščne klasifikatorje po intervalu bootstrap po klonih presegel le v enem od šestih zagonov. Tudi hipoteze H3 nismo potrdili, ker dodajanje fluksov k RNA in TCR napovedi ni zanesljivo izboljšalo pri nobenem od bolnikov.

Naknadna analiza je pokazala, da je pri oznaki ekspanzije iz surovih števil celic že velikost klona v krvi pred zdravljenjem napoved, ki je noben posamezen zagon modela ni presegel. Glavno spoznanje dela je zato, da AUC na vseh celicah pri tej nalogi ne meri, koliko informacije o ekspanziji nosita transkriptom ali metabolno stanje celice. Vhodne modalitete je treba primerjati ob izhodišču iz velikosti klona in med kloni, ki jih velikost klona ne loči. Pri takšni primerjavi je transkriptom dodatno informacijo nakazal pri klonih z eno celico pri enem bolniku, fluksi scFEA pa je niso pokazali.

Delo prispeva posodobljeno in vektorizirano izvedbo scFEA, s katero smo flukse ocenili za skoraj 70 000 celic, ter postopek vrednotenja, ki napoved modela TRIM primerja z velikostjo klona in negotovost oceni s ponovnim vzorčenjem po klonih. Največji težavi sta bili računska zahtevnost, zaradi katere smo model ovrednotili le pri dveh bolnikih, in nestabilnost učenja med semeni. Naslednji korak je preveriti, ali informacijo o ekspanziji pri klonih z eno celico nosijo metabolni geni, in model ovrednotiti pri več bolnikih z veliko ekspandiranimi kloni.

<!-- Izjava o avtorstvu (Vzorec 2026) že vsebuje alinejo o UI. Ali je potrebna še ločena izjava in kje,
     se dogovoriš z mentorjem; do takrat ostane tukaj kot neoštevilčen del, ki ni v kazalu. -->

::: {custom-style="Naslov brez stevilke"}
Izjava o uporabi umetne inteligence
:::


Pri pripravi dela sem uporabil orodje Claude (Anthropic) za pomoč pri strukturiranju besedila, jezikovnem pregledu in pisanju programske kode za obdelavo in analizo podatkov. Vse vsebine, vire in rezultate sem preveril in zanje v celoti odgovarjam.

# Literatura

::: {#refs}
:::
