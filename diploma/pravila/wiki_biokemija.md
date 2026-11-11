# Wiki FKKT: Navodila za pisanje zaključnih del (Biokemija)

Vir: <https://wiki.fkkt.uni-lj.si/index.php/Navodila_za_pisanje_zaklju%C4%8Dnih_del>
Verzija 1.2, posodobljena 3. 7. 2025 (avtor na strani ni naveden). Prebrano 11. 11. 2026.

**Za koga velja:** dodiplomski in magistrski študij **Biokemije**, za dela, oddana mentorju od 20. 6. 2022.
Kjer se razlikuje od splošnih pravil UL FKKT, **imajo ta navodila prednost**. Splošna navodila FKKT (2024) in Word vzorec veljata za vse ostalo (oblika strani, pisava, naslovnica, citiranje ...).

## Povezave na strani

| Dokument | Povezava |
|---|---|
| Splošna navodila FKKT za dipl./mag. delo (PDF, 2024) | povezava na wikiju vrne 404; **velja nova verzija 2026**, glej [README.md](README.md) |
| Vzorec diploma/magisterij (DOCX, 2024) | še deluje, a **velja Vzorec 2026**, glej [README.md](README.md) |
| Obrazec: Seznam nebesedilnih sestavin (XLSX, 2022c) | <http://web.fkkt.uni-lj.si/biokemija/Seznam%20nebesedilnih%20sestavin%202022c.xlsx> (11. 11. 2026 strežnik ni bil dosegljiv; po potrebi vprašaj mentorja ali tajništvo katedre) |

LaTeX predloge ni. O citatnem slogu, ključnih besedah, angleškem povzetku in UI stran ne pove ničesar (velja splošno, glej [README.md](README.md)).

## Pravila

**Poglavja (točna imena):** Uvod · Namen dela in hipoteze · Materiali in metode · Rezultati · Razprava · Zaključek · Literatura · (Priloge, neobvezno).
Rezultati in Razprava sta lahko združena, **če se tako dogovoriš z mentorjem**.

**Obseg:** diploma do **30 oštevilčenih strani brez literature in prilog**, dovoljeno odstopanje ±20 % (24–36 strani).
**Uvod** največ **1/3** dela (merjeno v oštevilčenih straneh brez literature).

**Povzetek:** ½ do 1 stran: področje, namen, metodološki pristop, ključni rezultati, kratek zaključek.

**Kratice:** slovenske kratice razložene slovensko (angleško le, če bi bil pomen sicer nejasen). V seznam **ne** gredo uveljavljene kratice (DNA, PCR, EDTA ...), kratice, uporabljene samo enkrat, in opisne razlage.

**Namen dela in hipoteze:**
- Pisano v **pretekliku** (kaj si želel ugotoviti, preden si začel delo).
- Lahko predvidene metode, **ne pa rezultatov**.
- Hipoteza = trditev, katere resničnosti pred poskusi ni bilo mogoče zanesljivo vedeti.
- Naslov po številu hipotez: »Namen dela in hipoteza« / »hipotezi« / »hipoteze«.

**Materiali in metode:** kratko. Reagenti abecedno s proizvajalcem in kataloško številko, oprema s proizvajalcem in modelom; brez držav proizvajalcev in znakov ®/™; rutinske metode brez podrobnega opisa; ključne metode podrobneje.

**Slike in tabele (nebesedilne sestavine):**
- Samo dve oznaki: **»slika«** (tudi sheme, fotografije, modeli) in **»tabela«**.
- Narediš jih **sam**. Že objavljene slike/tabele so dovoljene **samo v Uvodu**, z dovoljenjem založnika, razen če licenca (npr. CC) to izrecno dovoljuje; pogoje preveriš sam.
- Vse **v slovenščini**, tudi oznake osi in legende, z malimi začetnicami po slovenskem pravopisu.
- Strukture in modele narišeš iz javno dostopnih strukturnih zbirk.
- Oddaš **ločen list A4 s seznamom vseh nebesedilnih sestavin** in kako je bilo upoštevano avtorsko pravo (obrazec XLSX zgoraj).

**Priloge:** ne opisi metod in ne seznami materialov; npr. dodatne slike, zaporedja.

**Po končanem delu:** mentorju oddaš laboratorijske dnevnike in označene spremne datoteke; mentor jih hrani vsaj 10 let.

## Primerjava z našim povzetkom (README.md, 10. 11. 2026)

| Tema | README (splošna FKKT) | Wiki (Biokemija) | Ukrep |
|---|---|---|---|
| Poglavja | Rezultati in razprava skupaj ali ločeno | Ločeno; skupaj le po dogovoru z mentorjem | Vprašaj mentorja; trenutno `04_rezultati.md` = združeno |
| Namen | »Namen dela (in hipoteze)« 1–2 str. | »Namen dela in hipoteze«, preteklik, brez rezultatov | Naslov in čas v `02_namen.md` |
| Zaključek | »Zaključek / Sklepi« | »Zaključek« | ujema se |
| Obseg | ~30 strani | do 30 oštevilčenih str. brez lit. in prilog, ±20 % | natančneje |
| Uvod | ni omejitve | največ 1/3 | naš plan 8–9 str. je OK (≤10) |
| Povzetek | 3–5 ključnih besed, SLO + ANG | ½–1 stran, določena vsebina | dopolni |
| Kratice | brez SI enot in simbolov elementov | še brez uveljavljenih (DNA, PCR) in enkratnih | dopolni |
| Slike iz drugih virov | samo s pisnim dovoljenjem | **samo v Uvodu** in z dovoljenjem/licenco | strožje |
| Jezik slik | ni omenjeno | vse slovensko, tudi osi in legende | **pomembno**: vse grafe iz `analize/` je treba prevesti |
| Oznake | »slika«, tabele | samo »slika« in »tabela« | ujema se |
| Seznam nebesedilnih sestavin | ni omenjeno | ločen A4 list (XLSX obrazec) | **novo** |
| Lab. dnevniki | ni omenjeno | oddaš mentorju | pri nas: repozitorij + `lab_diary.md` |
| Povezava na Navodila 2024 | `.../Studij/pravilniki_in_obrazci/...` | `.../Studentski_referat/Pravilniki_studij/...` | skripta poskusi oba |
| Citiranje, oblika strani, UI | podrobno | ni omenjeno | velja README |

**Za računalniško diplomo** (brez mokrega laboratorija) smiselno prenesemo pravila za Materiale in metode: podatki (vir, GEO številka, število celic/bolnikov), programska oprema in paketi z **verzijami**, strojna oprema (Colab A100), nastavitve in semena; rutinske korake (npr. normalizacija scRNA) samo s sklicem.

## Kontrolni seznam (wiki FKKT, Biokemija)

Vir: <https://wiki.fkkt.uni-lj.si/index.php?title=Diplome_-_kontrolni_seznam> (prebrano 11. 11. 2026). Na vsako točko mora biti odgovor »da« ali utemeljitev za komisijo.

**Splošno**
- Številke poglavij brez pike na koncu (»1 Uvod«, »2.1 Podatki«).
- Tabela: naslov in kratek opis **nad** njo; slika: naslov in opis **pod** njo, na isti strani. Vsaka slika, tabela, enačba je omenjena v besedilu.
- Decimalna vejica; presledek med številom in enoto (tudi `25 °C`, `5 %`); nad 9999 ločilo tisočic (`146 000`).
- Številčni citati po vrstnem redu, znotraj stavka pred piko.
- Brez trpnika in tujk; žargon zamenjaj s slovenskimi izrazi. Ne »s pomočjo«, ampak »z«.
- Vsako novo poglavje se začne na desni (lihi) strani.
- Do 30 oštevilčenih strani (+20 %), brez literature in prilog. Priložen seznam nebesedilnih sestavin.

**Naslovnica:** po uradnem vzorcu 2024. **Povzetek in Abstract:** vsak največ 1 stran, angleški z  prevedenim naslovom, 3–5 ključnih besed.
**Kratice:** abecedno, razlaga slovensko, brez enot in kemijskih formul, vse uporabljene kratice.
**Uvod:** ≤1/3, vsaka trditev podprta z virom.
**Namen in hipoteze:** namen kot opis problema (ne povzetek), vsaj ena hipoteza kot trditev.
**Rezultati:** vsak prikazan rezultat opisan v besedilu s sklicem; »slika 3«, »tabela 2« z malo začetnico; **odstrani naslove in okvirje grafov, ki jih doda program** (naslov gre pod sliko).
**Razprava:** primerjava z drugimi raziskavami, **vsaj 5 objavljenih del**.
**Zaključek:** kratek in informativen, ne prepisan povzetek.
**Literatura:** dosledno; vsi avtorji (pri >10 prvi trije + »in sod.«); naslovi z malimi začetnicami (razen prve).
**Drugo:** kratica desno v zloženki (»celice HEK293«, »celice T CD8+«); zahvala (neobvezna) med izjavo in kazalom; seznami slik/tabel niso potrebni; v tabelah velika začetnica samo v glavi ali celih stavkih, števila poravnana na decimalno vejico; oznake slik slovensko, male začetnice.

Opomba: starejša stran wiki »Diplomsko delo« (2020) navaja povzetek 200–300 besed in ~50 000 znakov; novejši kontrolni seznam in navodila v1.2 imajo prednost.

## Primeri diplom

| Delo | Zakaj je uporabno |
|---|---|
| F. Petrovič: *Modeliranje metabolizma duktalnega adenokarcinoma na ravni genoma*, diplomsko delo, UL FKKT, 2025, mentor M. Moškon. [RUL](https://repozitorij.uni-lj.si/IzpisGradiva.php?id=171844) | Isti mentor, FKKT, ftINIT in csGEM iz transkriptomskih podatkov: najboljši vzor za strukturo in obliko |
| D. Vaupotič: *Modeliranje celične rasti*, diplomsko delo, UL FKKT, 2019. [RUL](https://repozitorij.uni-lj.si/IzpisGradiva.php?id=110310) | Računalniška diploma na FKKT (brez laboratorija) |
| J. Nastran: *Integration and analysis of scRNAseq data using genome-scale metabolic models*, magistrsko delo, mentor M. Moškon. [RUL](https://repozitorij.uni-lj.si/IzpisGradiva.php?id=155937) | Vsebinsko zelo blizu (scRNA + GEM), za literaturo in metode; fakulteta ni preverjena, oblika verjetno ni po FKKT |

## Petrovič 2025: kaj pokaže dejanska diploma pri Moškonu

F. Petrovič, *Modeliranje metabolizma duktalnega adenokarcinoma na ravni genoma*, UL FKKT, Biokemija, mentor prof. dr. M. Moškon, 2025 (67 strani PDF, prebrano 11. 11. 2026). Oddana je bila pred pravili 2026.

| Del | Kako je narejeno | Za nas |
|---|---|---|
| Naslovnica, notranja naslovnica, izjava | po Vzorcu (izjava še brez alineje o UI) | uporabimo Vzorec 2026 |
| Prazne strani | vsak del prvega odseka in vsako poglavje se začne na lihi strani | naša gradnja to naredi samodejno |
| Glava | »F. Petrovič, *naslov*« na vsaki strani drugega odseka | neobvezno; dodamo, če želiš |
| Povzetek / Abstract | vsak na svoji strani, naslov dela nad njim, 4 ključne besede | enako |
| Kratice | slovenska razlaga + »(angl. ...)«; vključene tudi DNA, ATP | po kontrolnem seznamu Biokemije DNA, RNA ipd. izpustimo |
| Uvod | str. 1–6 (pod 1/3) | ciljamo 8–9 strani |
| Namen | »Namen dela«, pol strani, **brez hipotez** | pravila Biokemije zahtevajo hipoteze, zato jih obdržimo |
| Materiali in metode | str. 9–17: podatki, model Human1, ftINIT, pragovi, statistični testi, vizualizacija | podoben vrstni red: podatki → predobdelava → modeli → evalvacija |
| Rezultati in razprava | **združeno** v eno poglavje, str. 19–31 | po pravilih Biokemije ločeno (Moškon je pri Petroviču združeno sprejel) |
| Zaključek | str. 33 | 1–2 strani |
| Literatura | str. 35–40, 53 virov; imena revij večinoma izpisana v celoti (ne po CAS) | držimo se Navodil: kratice revij in DOI |
| Priloge | dodatni grafi (P1, P2 ...) | po potrebi |
| Slike | »Slika 7: ...« pod sliko, tabela z naslovom nad njo | enako |

Petrovičeva diploma je dober vzor za **obseg in strukturo računalniške diplome pri Moškonu**, ni pa vir pravil: kjer se razlikuje od Navodil 2026 ali pravil Biokemije, sledimo pravilom.
