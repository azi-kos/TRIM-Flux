# Orodja za delo z literaturo (poceni in zanesljivo)

Stanje 10. 11. 2026. Cilj: pokriti veliko znanstvene literature (scRNA-seq, TCR, imunoterapija HNSCC, VAE/multimodalno učenje, metabolni modeli GEM/FBA), ne da bi Claude bral cele članke po nepotrebnem, in da je vsak citat preverljiv.

## Kratek seznam (priporočeno)

| Prioriteta | Orodje | Za kaj | Cena | Stanje pri tebi |
|---|---|---|---|---|
| 1 | **PubMed** connector | iskanje, povzetki, metapodatki, polna besedila iz PMC (odprti dostop), sorodni članki, preverjanje citatov | brezplačno, brez računa | ni dodan |
| 1 | **Zotero** (namizni) + Zotero Connector (brskalnik) + **Better BibTeX** | ena knjižnica virov, avtomatski `.bib`, citati v Wordu ali Pandoc, slog `fkkt-ul.csl` | brezplačno | namesti sam |
| 1 | **Google Drive** connector | sem odložiš PDF-je, ki jih dobiš prek UL dostopa (DiKUL), da jih lahko preberem | že povezan | ✓ |
| 2 | **bioRxiv/medRxiv** connector | najnovejši preprinti (scFEA, METAFlux, TRIM-podobne metode) | brezplačno | ni dodan |
| 2 | **alphaXiv** connector | arXiv: VAE, multimodalno učenje, ML metode | brezplačen račun | ni dodan |
| 3 | Consensus | hitri »kaj pravi literatura o X« odgovori s citati | omejen brezplačni paket | ni dodan, neobvezno |
| — | Scite, Elicit, Wiley Scholar Gateway, SciSpace | dobra, a večinoma plačljiva ali podvajajo PubMed | plačljivo / delno | **ne priporočam** za diplomo |

Zotero connectorja v imeniku Claude ni. Zotero zato uporabljamo kot vir resnice za reference: iz njega izvoziš `literatura.bib` (Better BibTeX, samodejni izvoz) v repo, jaz pa berem in citiram samo iz te datoteke.

DiKUL (digitalna knjižnica UL, prijava z UL računom) ti da dostop do večine plačljivih revij. Tega dostopa jaz nimam, zato zaprte PDF-je preneseš ti.

Za diplomo ne uporabljamo »deep research« ali več vzporednih agentov: so dragi, za ~40–60 virov pa niso potrebni.

## Postopek: najprej povzetki, nato polna besedila

1. **Seznam vprašanj po poglavjih.** Za vsako poglavje 3–6 trditev, ki jih moramo podpreti z virom (npr. »klonska ekspanzija T-celic napoveduje odziv na anti-PD-1 pri HNSCC«).
2. **Iskanje samo po povzetkih (PubMed, bioRxiv, alphaXiv).** Za vsako trditev 3–5 kandidatov; berem samo naslov, povzetek in metapodatke. To je približno 10-krat cenejše od branja celih člankov.
3. **Razvrstitev v tri skupine:**
   - **A, jedro (~10–15 virov):** TRIM, podatki GSE200996, scFEA, METAFlux, Human-GEM, ftINIT/RAVEN, scVI/VAE, ključni članki o TCR in imunoterapiji. Te preberemo v celoti: odprti dostop prek PMC, ostalo kot PDF na Drive.
   - **B, podporni (~20–30):** dovolj sta povzetek in ena preverjena trditev.
   - **C, ozadje (~10):** citiramo pregledni članek namesto desetih primarnih.
4. **Ena kartica na vir** v `literatura/opombe/<citekey>.md`: 5–10 vrstic (kaj trdi, ključne številke, kje v diplomi ga citiramo, ali je preverjeno v originalu). Ko je kartica napisana, članka ni več treba brati: pri pisanju uporabljam kartice, ne PDF-jev.
5. **Preverjanje pred oddajo:** vsak `[n]` v besedilu se poveže z DOI v `literatura.bib`. Z orodjem PubMed (`lookup_article_by_citation`) ali prek DOI preverimo, da vir obstaja in pravi to, kar trdimo. Smernice UL izrecno prepovedujejo izmišljene vire.

## Kje in kako pišemo

Diplomo pišemo **v repozitoriju kot Markdown** (`diploma/poglavja/`, ena datoteka na poglavje), citati so oblike `[@citekey]`, in jo s Pandoc pretvorimo v Word (`bash diploma/zgradi.sh`). Postopek je opisan v [`diploma/README.md`](../README.md).

Prednosti: lahko pišem in popravljam neposredno v repozitoriju, vsaka sprememba je vidna v gitu, končni izgled pa da predloga FKKT. Zadnji korak (naslovnica, izjava, PDF/A) narediš v Wordu.

## Kako dodaš priporočene connectorje

V claude.ai: Nastavitve → Connectors → Browse → poišči **PubMed**, **bioRxiv**, **alphaXiv** (po želji Consensus) → Connect. PubMed in bioRxiv ne potrebujeta računa.
