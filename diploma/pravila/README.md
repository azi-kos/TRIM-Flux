# Pravila in predloge za diplomsko delo (UL FKKT)

Zbrano 10. 11. 2026, **posodobljeno 11. 11. 2026 na pravila, ki veljajo od 1. 10. 2026** (originali prebrani neposredno s strani FKKT). Velja za diplomsko delo 1. stopnje na UL FKKT. **Za študij Biokemije** veljajo še strožja navodila z wiki FKKT, ki imajo prednost pred splošnimi: glej [`wiki_biokemija.md`](wiki_biokemija.md). Ko se kaj ne ujema, ima prednost dogovor z mentorjem (prof. Moškon) in najnovejša verzija na spletni strani FKKT.

> **Kaj se je spremenilo 25. 9. 2026 (velja od 1. 10. 2026):** tema velja **2 leti** (prej 10 mesecev); končni PDF/A v Studis **vsaj 5 delovnih dni** pred zagovorom (prej 7); potrdilo knjižnice je treba pridobiti pred nalaganjem končne verzije; mentor ima lahko največ 3 študente; **Izjava o avtorstvu ima novo alinejo o UI**: »sem v primeru uporabe orodij umetne inteligence le-ta uporabljal/-a zgolj kot pomoč pri oblikovanju besedila ter zbiranju in obdelavi podatkov«; nova Navodila in Vzorec 2026.

**Originalne datoteke** (niso v repozitoriju) prenese skripta:

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
| `wiki_biokemija.md` | navodila z wiki FKKT za študij Biokemije (v1.2, 3. 7. 2025) in primerjava s tem povzetkom |

## 1. Uradni dokumenti

| # | Dokument | Verzija | Povezava |
|---|---|---|---|
| 1 | **Pravilnik o diplomskem delu 1. stopnje, prečiščeno besedilo** | 25. 9. 2026, velja od 1. 10. 2026 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/habilitacije/PRAVILNIK_O_DIPLOMSKEM_DELU_PRVOSTOPENJSKIH_STUDIJSKIH_PROGRAMOV_FAKULTETE_ZA_KEMIJO_IN_KEMIJSKO_TEHNOLOGIJO_UNIVERZE_V_LJUBLJANI.pdf) |
| 1a | Pravilnik o spremembah (kaj je novo) | 25. 9. 2026 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/o_fakulteti/Pravilnik_o_spremembah_Pravilnika_o_diplomskem_delu_prvostopenjskih_studijskih_programov_Fakultete_za_kemijo_in_kemijsko_tehnologijo_Univerze_v_Ljubl.pdf) |
| 2 | **Navodila za izdelavo diplomskega ali magistrskega dela** (Priloga 3) | 2026, čistopis | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/habilitacije/Navodila_za_izdelavo_dipl_mag_dela_2026_-_cistopis.docx) |
| 3 | **Vzorec diplomskega dela (Word, Priloga 2)** | 2026, čistopis | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/habilitacije/Vzorec_diploma_magisterij_2026_-_cistopis.docx) |
| 4 | Template (angleška Word predloga) | 2026 | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Template_diploma_master_2026.docx) |
| 5 | Rokovnik dejavnosti za diplomsko delo | posodobljen 13. 10. 2026 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/ROKOVNIK_DEJAVNOSTI_ZA_DIPLOMSKO_DELO.pdf) |
| 6 | Dispozicija diplomskega dela (samo za delo izven UL FKKT) | 2025 | [DOCX](https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Dispozicija_diplomskega_dela_studenta.docx) |

Vse povezave so s strani [Univerzitetni študijski programi](https://fkkt.uni-lj.si/studij-in-studentje/univerzitetni-studijski-programi) (preverjeno 11. 11. 2026). Starejše povezave (Pravilnik 2023, Navodila in Vzorec 2024 na starih poteh) vrnejo 404 in **ne veljajo več**.
| 7 | UL in UI: izjava in smernice (rektor) | 19. 9. 2023 | [PDF](https://fkkt.uni-lj.si/fileadmin/user_upload/UL_in_UI_izjava_in_smernice_za_objavo_1_.pdf) |
| 8 | Smernice UL za uporabo UI (novejše, celotna UL) | junij 2026 | [PDF (kopija na TEOF)](https://www.teof.uni-lj.si/uploads/UL26_smernice-UI(1).pdf) |
| 9 | Odgovorno in etično raziskovanje (UL) | | [PDF](https://www.uni-lj.si/assets/Sluzba-za-raziskovalno-dejavnost/Etika-in-integriteta/Odgovorno-in-eticno-raziskovanje.pdf) |
| 10 | Statut UL, 133. člen (oddaja, repozitorij RUL, preverjanje podobnosti) | Ur. l. RS 49/2022 | [Uradni list](https://www.uradni-list.si/glasilo-uradni-list-rs/vsebina/2022-01-1002) |

Uradne LaTeX predloge FKKT **ni**. Obstajajo le predloge drugih članic UL (npr. FE, FRI) na [Overleaf (uni-lj)](https://www.overleaf.com/latex/templates/tagged/uni-lj), ki pa ne sledijo oblikovanju FKKT. Če pišemo v LaTeX/Markdown, moramo oblikovanje FKKT prenesti sami (glej točko 3).

Starejša navodila (`Navodila_izdelava_diploma_magisterij_koncna.pdf`, verzija 2018) zahtevajo tiskane izvode v trdih platnicah in lastnoročni podpis. **To ne velja več**: od 2024 se oddaja samo elektronsko (PDF/A, digitalno podpisana izjava).

## 2. Postopek in roki (Pravilnik 2026 + Rokovnik)

1. **Mentor** se dodeli prek ankete ob vpisu v 3. letnik (zavezujoče), seznami do 30. 11.
2. **Diplomski tečaj** v spletni učilnici (pregled vsebin + vprašalnik) je pogoj za odobritev teme.
3. **Prijava teme** v Studisu (obrazec »Prošnja za določitev teme diplomskega dela«) **do 1. marca**. Mentor potrdi in predlaga člane komisije, prodekan odobri (prvi teden marca).
4. **Tema velja 2 leti** od odobritve (18. člen: delo je treba predložiti mentorju najkasneje 2 leti po odobritvi teme; Rokovnik: v tem času ga pripraviti in zagovarjati). Sicer tema zapade.
5. **Naslov** se lahko spremeni do oddaje komisiji (mentorjevo soglasje, e-pošta referatu).
6. **Pregled:** najprej mentor, nato člani komisije; komisija vrne pripombe v **21 dneh**.
7. **Pred nalaganjem končne verzije potrdilo knjižnice** o poravnanih obveznostih. **Končna verzija (PDF/A) v Studis vsaj 5 delovnih dni pred zagovorom** (nalaganje omogoči mentor po odobritvi), z digitalno podpisano Izjavo o avtorstvu. Sledi **preverjanje podobnosti**; mentor in komisija potrdita poročilo, sicer popravek in ponovno nalaganje.
8. **Prijava zagovora (datum in ura, usklajena s komisijo) + oddaja potrdila knjižnice v referat vsaj 3 delovne dni pred zagovorom.**
9. **Zagovor:** komisija = mentor + vsaj en član z učiteljskim nazivom UL; somentor sodeluje, ne ocenjuje. Ocena upošteva ustno predstavitev, kakovost **prve oddane verzije** in upoštevanje popravkov, poznavanje področja in teorije. Z odliko: povprečje ≥ 9,0 in diploma 10.
10. **Objava:** elektronska verzija gre v Repozitorij UL (RUL); s podpisom študent podeli UL neizključno pravico javne objave. Začasna nedostopnost (npr. zaradi objave članka) je mogoča po pravilniku UL, zaprositi je treba pred oddajo.

> Ker se ocenjuje tudi **prva oddana verzija**, naj bo prva verzija, ki gre mentorju, že jezikovno in oblikovno urejena.

## 3. Oblika (Navodila 2026 + Vzorec 2026)

- **A4, robovi 3 cm** povsod; **Times New Roman ali Arial 12, razmik 1,15**.
- **Obseg:** okvirno **~30 strani** (diploma), magisterij 50–60.
- **Vzorec 2026 (dejanski slogi):** besedilo TNR 12, obojestransko poravnano, 12 pt za odstavkom; naslov 1. ravni 16 krepko (30 pt za njim), 2. raven 14 krepko (24 pt pred), 3. raven 12 krepko ležeče; napisi slik sredinsko; **vsako poglavje se začne na lihi strani** (prelomi razdelkov »liha stran«). Naša gradnja (`zgradi.sh`) to posnema.
- **Jezik:** slovenščina (angleščina samo na odobreno prošnjo + razširjen slovenski povzetek 3–5 strani). Isti pojem vedno enako, brez nepotrebnih tujk; nov strokovni izraz brez slovenskega ustreznika: ob prvi pojavitvi angleški izraz v oklepaju, npr. »klonska ekspanzija (angl. *clonal expansion*)«.
- **Naslovi:** največ 3 ravni. 1. raven 16 krepko, 2. raven 14 krepko, 3. raven 12 krepko-ležeče.
- **Oštevilčenje:** prvi odsek brez številk ali male rimske; od Uvoda naprej arabske, spodaj na sredini. Glava neobvezna.

### Vrstni red

Prvi odsek:
1. Naslovnica: UNIVERZA V LJUBLJANI / FAKULTETA ZA KEMIJO IN KEMIJSKO TEHNOLOGIJO / DIPLOMSKO DELO (20 pt) / ime in priimek / Ljubljana, leto. Sredinsko, 16 pt.
2. Notranja naslovnica: UNIVERZA V LJUBLJANI / FAKULTETA ZA KEMIJO IN KEMIJSKO TEHNOLOGIJO / UNIVERZITETNI ŠTUDIJSKI PROGRAM BIOKEMIJA / **naslov dela (krepko)** / DIPLOMSKO DELO (20 pt) / ime in priimek / MENTOR: prof. dr. Miha Moškon / Ljubljana, leto.
3. Izjava o avtorstvu (naslov 18 krepko, »diplomskega dela« 16), elektronsko podpisana; alineje: lastno delo pod mentorstvom, citiranje po navodilih, plagiatorstvo kaznivo (ZASP-UPB3), slovnična in oblikovna korektnost, afiliacija UL FKKT ob objavi, **UI samo kot pomoč pri oblikovanju besedila ter zbiranju in obdelavi podatkov** (novo 2026).
4. Zahvala (neobvezna, ležeče).
5. **Povzetek + ključne besede (3–5)** v slovenščini in angleščini (»Abstract«, »Keywords:«), z naslovom dela; brez slik, razumljiv samostojno.
6. Kazalo (naslov »Kazalo« 16 krepko; 3 ravni; neoštevilčeni deli, kot so povzetek in izjava, ne sodijo vanj).
7. »Seznam uporabljenih kratic in simbolov« (16 krepko; abecedno ali tematsko; kratice razložene slovensko in angleško, kjer je potrebno; brez SI enot, simbolov elementov in izpeljank iz drugih kratic).

Drugi odsek (poglavja se lahko prilagodijo v dogovoru z mentorjem):
1. **Uvod** (pregled literature, problem, cilji)
2. **Namen dela** (in **hipoteze**), 1–2 strani
3. **Eksperimentalni del / Materiali in metode**: dovolj natančno za ponovitev; logični, ne kronološki vrstni red; vsaka metoda iz literature z referenco; za nas: podatki, programska oprema, verzije, strojna oprema, nastavitve.
4. **Rezultati in razprava** (skupaj ali ločeno): primerjava s hipotezami in literaturo.
5. **Zaključek / Sklepi**: ocena hipotez, težave, nadaljnje delo.
6. **Literatura** (oštevilčeno poglavje, npr. »6 Literatura«), nato po potrebi **Priloge** (»7 Priloge«).

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

**Nova alineja v Izjavi o avtorstvu (od 1. 10. 2026)** – to podpišeš z digitalnim podpisom:
> »sem v primeru uporabe orodij umetne inteligence le-ta uporabljal/-a zgolj kot pomoč pri oblikovanju besedila ter zbiranju in obdelavi podatkov.«

To je strožje od smernic UL: UI sme pomagati pri **oblikovanju** besedila (jezik, struktura, preglednost) ter pri **zbiranju in obdelavi podatkov** (iskanje literature, koda za analizo). Vsebine (argumentov, interpretacije, zaključkov) ne sme pisati namesto tebe.

**Za nas to pomeni:**
1. Uporabo UI za pisanje diplome **vnaprej uskladi z mentorjem** (smernice 2026 to izrecno priporočajo za zaključna dela) in se dogovori, kam gre izjava.
2. Delovni tok, skladen z izjavo: ti napišeš vsebino (lahko v alinejah ali grobo); Claude pomaga pri strukturi, jeziku in obliki, poišče in preveri vire ter piše kodo za analizo in slike. Interpretacijo in zaključke oblikuješ ti. Vsak citat preverimo v originalnem viru (DOI).
3. Na koncu dodamo kratko **izjavo o uporabi UI**, npr.: »Pri pripravi dela sem uporabil orodje Claude (Anthropic, november 2026–…) za pomoč pri strukturiranju besedila, jezikovnem pregledu in pisanju programske kode za analizo. Vse vsebine, vire in rezultate sem preveril in zanje odgovarjam.«
