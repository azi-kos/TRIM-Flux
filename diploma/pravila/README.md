# Pravila in predloge za diplomsko delo (UL FKKT)

Zbrano 10. 11. 2026. Velja za diplomsko delo 1. stopnje na UL FKKT. Ko se kaj ne ujema, ima prednost dogovor z mentorjem (prof. Moškon) in najnovejša verzija na spletni strani FKKT.

**Originalne datoteke:** v tem okolju strežnik fkkt.uni-lj.si ni dosegljiv, zato originalov ni v repozitoriju. Prenese jih skripta:

```bash
bash diploma/pravila/prenesi_dokumente.sh   # -> diploma/pravila/originali/
```

Datoteke v mapi:

| Datoteka | Kaj je |
|---|---|
| `README.md` | ta povzetek pravil + povezave na uradne dokumente |
| `prenesi_dokumente.sh` | prenos vseh originalov (PDF, DOCX) |
| `fkkt-ul.csl` | neuradni citatni slog za Zotero/Pandoc, ki točno sledi primeru iz navodil FKKT |
| `orodja_za_literaturo.md` | priporočeni connectorji, orodja in poceni postopek dela z literaturo |

## 1. Uradni dokumenti

| # | Dokument | Verzija | Povezava |
|---|---|---|---|
| 1 | Pravilnik o diplomskem delu (čistopis) | 20. 10. 2023, velja od 1. 1. 2024 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/cistopis_PRAVILNIK_O_DIPLOMSKEM_DELU_okt_2023.pdf) |
| 2 | Navodila za izdelavo diplomskega in magistrskega dela | 2024 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki_in_obrazci/Navodila_za_izdelavo_dipl_mag_dela_2024.pdf) |
| 3 | **Vzorec (Word predloga)** diploma/magisterij | 2024 | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Vzorec_diploma_magisterij_2024.docx) |
| 4 | Template (angleška Word predloga) | 2024 | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Template_diploma_master_2024.docx) |
| 5 | Rokovnik dejavnosti za diplomsko delo | aktualen | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/ROKOVNIK_DEJAVNOSTI_ZA_DIPLOMSKO_DELO.pdf) |
| 6 | Pravilnik 1. stopnja | aktualen | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/Pravilnik_1_stopnja.pdf) |
| 7 | UL in UI: izjava in smernice (rektor) | 19. 9. 2023 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/UL_in_UI_izjava_in_smernice_za_objavo_1_.pdf) |
| 8 | Smernice UL za uporabo UI (novejše, celotna UL) | junij 2026 | [PDF (kopija na TEOF)](https://www.teof.uni-lj.si/uploads/UL26_smernice-UI(1).pdf) |
| 9 | Odgovorno in etično raziskovanje (UL) | | [PDF](https://www.uni-lj.si/assets/Sluzba-za-raziskovalno-dejavnost/Etika-in-integriteta/Odgovorno-in-eticno-raziskovanje.pdf) |
| 10 | Statut UL, 133. člen (oddaja, repozitorij RUL, preverjanje podobnosti) | Ur. l. RS 49/2022 | [Uradni list](https://www.uradni-list.si/glasilo-uradni-list-rs/vsebina/2022-01-1002) |

Uradne LaTeX predloge FKKT **ni**. Obstajajo le predloge drugih članic UL (npr. FE, FRI) na [Overleaf (uni-lj)](https://www.overleaf.com/latex/templates/tagged/uni-lj), ki pa ne sledijo oblikovanju FKKT. Če pišemo v LaTeX/Markdown, moramo oblikovanje FKKT prenesti sami (glej točko 3).

Starejša navodila (`Navodila_izdelava_diploma_magisterij_koncna.pdf`, verzija 2018) zahtevajo tiskane izvode v trdih platnicah in lastnoročni podpis. **To ne velja več**: od 2024 se oddaja samo elektronsko (PDF/A, digitalno podpisana izjava).

## 2. Postopek in roki (Pravilnik 2023 + Rokovnik)

1. **Mentor** se dodeli prek ankete ob vpisu v 3. letnik (zavezujoče), seznami do 30. 11.
2. **Diplomski tečaj** v spletni učilnici (pregled vsebin + vprašalnik) je pogoj za odobritev teme.
3. **Prijava teme** v Studisu (obrazec »Prošnja za določitev teme diplomskega dela«) **do 1. marca**. Mentor potrdi in predlaga člane komisije, prodekan odobri (prvi teden marca).
4. **Tema velja 10 mesecev** od odobritve. Podaljšanje: največ 1 mesec na utemeljeno prošnjo; če je prva verzija oddana pravočasno, z mentorjevim soglasjem do konca študijskega leta. Sicer tema zapade.
5. **Naslov** se lahko spremeni do oddaje komisiji (mentorjevo soglasje, e-pošta referatu).
6. **Pregled:** najprej mentor, nato člani komisije; komisija vrne pripombe v **21 dneh**.
7. **Končna verzija (PDF/A) v Studis vsaj 7 delovnih dni pred zagovorom**, z digitalno podpisano Izjavo o avtorstvu. Sledi **preverjanje podobnosti**; mentor in komisija potrdita poročilo, sicer popravek in ponovno nalaganje.
8. **Prijava zagovora + potrdilo knjižnice vsaj 3 delovne dni pred zagovorom.**
9. **Zagovor:** komisija = mentor + vsaj en član z učiteljskim nazivom UL; somentor sodeluje, ne ocenjuje. Ocena upošteva ustno predstavitev, kakovost **prve oddane verzije** in upoštevanje popravkov, poznavanje področja in teorije. Z odliko: povprečje ≥ 9,0 in diploma 10.
10. **Objava:** elektronska verzija gre v Repozitorij UL (RUL); s podpisom študent podeli UL neizključno pravico javne objave. Začasna nedostopnost (npr. zaradi objave članka) je mogoča po pravilniku UL, zaprositi je treba pred oddajo.

> Ker se ocenjuje tudi **prva oddana verzija**, naj bo prva verzija, ki gre mentorju, že jezikovno in oblikovno urejena.

## 3. Oblika (Navodila 2024)

- **A4, robovi 3 cm** povsod; **Times New Roman ali Arial 12, razmik 1,15**.
- **Obseg:** okvirno **~30 strani** (diploma), magisterij 50–60.
- **Jezik:** slovenščina (angleščina samo na odobreno prošnjo + razširjen slovenski povzetek 3–5 strani). Isti pojem vedno enako, brez nepotrebnih tujk; nov strokovni izraz brez slovenskega ustreznika: ob prvi pojavitvi angleški izraz v oklepaju, npr. »klonska ekspanzija (angl. *clonal expansion*)«.
- **Naslovi:** največ 3 ravni. 1. raven 16 krepko, 2. raven 14 krepko, 3. raven 12 krepko-ležeče.
- **Oštevilčenje:** prvi odsek brez številk ali male rimske; od Uvoda naprej arabske, spodaj na sredini. Glava neobvezna.

### Vrstni red

Prvi odsek:
1. Naslovnica: UNIVERZA V LJUBLJANI / FAKULTETA ZA KEMIJO IN KEMIJSKO TEHNOLOGIJO / ime in priimek / DIPLOMSKO DELO (20 pt) / Ljubljana, leto. Vse velike črke, sredinsko, 16 pt.
2. Notranja naslovnica: še študijski program, **naslov dela (krepko)**, mentor (in somentor).
3. Izjava o avtorstvu (naslov 18 krepko), digitalno podpisana; vključuje zavedanje o kaznivosti plagiatorstva (ZASP) in afiliacijo UL FKKT ob objavi.
4. Zahvala (neobvezna, ležeče).
5. **Povzetek + ključne besede (3–5)** v slovenščini in angleščini (»Abstract«, »Keywords:«), z naslovom dela; brez slik, razumljiv samostojno.
6. Kazalo (3 ravni).
7. Seznam kratic in simbolov (abecedno ali tematsko; brez SI enot in simbolov elementov).

Drugi odsek (poglavja se lahko prilagodijo v dogovoru z mentorjem):
1. **Uvod** (pregled literature, problem, cilji)
2. **Namen dela** (in **hipoteze**), 1–2 strani
3. **Eksperimentalni del / Materiali in metode**: dovolj natančno za ponovitev; logični, ne kronološki vrstni red; vsaka metoda iz literature z referenco; za nas: podatki, programska oprema, verzije, strojna oprema, nastavitve.
4. **Rezultati in razprava** (skupaj ali ločeno): primerjava s hipotezami in literaturo.
5. **Zaključek / Sklepi**: ocena hipotez, težave, nadaljnje delo.
6. **Literatura**.

### Slike, tabele, enačbe

- Vse oštevilčeno (enotno, npr. »slika 7« ali »slika 2.1«); vsako omeniti v besedilu na prvem smiselnem mestu.
- **Slik iz drugih virov ni dovoljeno kopirati brez pisnega dovoljenja** imetnika avtorskih pravic. Priporočilo: vse slike narišemo sami (arhitektura TRIM, shema poteka itd.) in navedemo »prirejeno po [n]«.

### Citiranje

- V besedilu **številke v oglatem oklepaju**, **pred ločilom**, zaporedno po prvi pojavitvi: »… klonsko ekspanzijo [3].«
- Seznam literature v istem vrstnem redu, z naslovi del. Kratice revij po CAS Source Index.
- Članek:
  `[1] M. P. C. Margues, P. Fernandes, J. M. S. Cabral, P. Žnidaršič Plazl, I. Plazl: Continuous steroid biotransformations in microchannel reactors. New Biotechnol. 2012, 29, 227–234. DOI: 10.1016/j.nbt.2011.10.001`
  (kratica revije in letnik ležeče)
- Knjiga: `Avtorji: Naslov. Kraj: Založba, leto, str. xx–xy.`
- Poglavje: `Avtorji: Naslov poglavja. V: Naslov knjige, izdaja, zvezek, urednik (ur.), Kraj: Založba, leto, str. xx–xy.`
- Spletne strani, patenti, diplomska dela: po primerih v Word vzorcu.
- Ta slog je implementiran v [`fkkt-ul.csl`](fkkt-ul.csl) (preverjeno, izpiše natanko zgornji primer).

## 4. Plagiatorstvo

- Plagiatorstvo (tuje misli/besedilo/slike kot lastne, prepisovanje brez vira) vodi v disciplinski postopek; po ZASP je kaznivo.
- Vsa dela gredo skozi **program za preverjanje podobnosti**; poročilo pregledata mentor in komisija.
- Tudi parafraza tujega besedila (ročna ali z UI) brez navedbe vira je plagiat.

## 5. Umetna inteligenca (UI)

Po smernicah UL (2023 in 2026):

**Dovoljeno** (če mentor ne določi drugače): jezikovni pregled lastnega besedila, prevod in parafraza **lastnega** besedila, orientacija po temi, iskalni izrazi, pomoč pri strukturiranju, programiranju, analizi podatkov in vizualizacijah.

**Prepovedano:**
- oddaja pretežno UI-generiranega besedila kot lastnega ali prepuščanje osrednjih argumentov, interpretacije in zaključkov UI;
- dobesedno kopiranje UI-besedila brez navedbe;
- parafraziranje tujih besedil z UI, da se izognemo citiranju;
- **izmišljeni viri, citati, podatki ali rezultati** (vsak vir mora biti preverjen v originalu);
- vnos zaupnih, osebnih ali avtorsko varovanih gradiv v javna UI orodja;
- navajanje UI kot (so)avtorja.

**Razkritje:** obvezno, kadar UI pomembno vpliva na vsebino, strukturo, analizo, rezultate ali zaključke, ali kadar to zahteva mentor. Navede se orodje (in verzija), datum, namen, deli dela in način preverjanja. Mesto: izjava o uporabi UI, metode, opomba ali priloga, po dogovoru z mentorjem. Detektorji UI niso samostojen dokaz kršitve, študent pa mora znati pojasniti svoj postopek.

**Za nas to pomeni:**
1. Uporabo UI za pisanje diplome **vnaprej uskladi z mentorjem** (smernice 2026 to izrecno priporočajo za zaključna dela) in se dogovori, kam gre izjava.
2. Delovni tok: Claude pripravi osnutke, strukturo in kodo; ti vsebino, interpretacijo in zaključke preveriš, preoblikuješ in jih oblikuješ po svoje. Vsak citat preverimo v originalnem viru (DOI).
3. Na koncu dodamo kratko **izjavo o uporabi UI**, npr.: »Pri pripravi dela sem uporabil orodje Claude (Anthropic, november 2026–…) za pomoč pri strukturiranju besedila, jezikovnem pregledu in pisanju programske kode za analizo. Vse vsebine, vire in rezultate sem preveril in zanje odgovarjam.«
