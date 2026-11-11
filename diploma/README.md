# Diplomsko delo – TRIM-Flux

Vse, kar je povezano s pisanjem diplome, je v tej mapi. Koda in rezultati projekta ostanejo nespremenjeni v preostalem repozitoriju.

```
diploma/
├── metadata.yaml          naslov, avtor, mentor (PREVERI podatke)
├── poglavja/              besedilo, ena Markdown datoteka na poglavje (vrstni red = ime datoteke)
├── literatura/
│   ├── literatura.bib     reference (izvoz iz Zotera, Better BibTeX)
│   ├── kandidati.md       seznam virov, ki jih še potrebujemo
│   └── opombe/            ena kratka kartica na vir (_predloga.md)
├── slike/                 lastne slike (PNG/SVG), vstavi: ![Slika 1: opis](slike/ime.png)
├── predloga/              Word referenca po Vzorcu FKKT 2026 + razdelki.py (naslovnici, izjava, lihe strani)
├── pravila/               pravila FKKT/UL, citatni slog fkkt-ul.csl, orodja za literaturo
└── zgradi.sh              Markdown -> izhod/diploma.docx
```

## Gradnja

```bash
bash diploma/zgradi.sh                 # potrebuje pandoc (>= 3) in python-docx
```

Izhod sledi **Vzorcu UL FKKT 2026**: naslovnica, notranja naslovnica in izjava o avtorstvu (z alinejo o UI) se naredijo iz `metadata.yaml`; vsak del in vsako poglavje se začne na lihi strani; prvi odsek je brez številk strani, od Uvoda naprej so številke spodaj na sredini. Word ob odprtju vpraša, ali naj posodobi polja: potrdi, da se izpolni kazalo. Na koncu izvoz v PDF/A in elektronski podpis izjave.

V `metadata.yaml` nastavi `spol: m` ali `spol: ž`, da izjava uporabi pravilne oblike (podpisani/podpisana ...).

## Pisanje

- Citat: `... klonsko ekspanzijo [@kljuc2023].` V Wordu postane `[1]`, pred piko, oštevilčeno po vrstnem redu.
- Komentarji `<!-- ... -->` so delovne opombe in se ne izpišejo.
- Obseg: do 30 oštevilčenih strani brez literature in prilog (±20 %). Okvirno: Uvod 8–9 (največ 1/3), Namen 1–2, Materiali in metode 7–8, Rezultati 5, Razprava 4, Zaključek 1–2.
- Vrstni red poglavij po navodilih za Biokemijo: Uvod, Namen dela in hipoteze, Materiali in metode, Rezultati, Razprava, Zaključek, Literatura.
- Pred oddajo preveri [kontrolni seznam](pravila/wiki_biokemija.md#kontrolni-seznam-wiki-fkkt-biokemija). V Wordu nastavi, da se vsako poglavje začne na lihi strani (prelom razdelka »Liha stran«), in priloži seznam nebesedilnih sestavin (XLSX).
- Številke v rezultatih vpisujemo iz shranjenih logov v `results/logs/`, ne na pamet.
- Pravila oblike, citiranja in uporabe UI so povzeta v [`pravila/README.md`](pravila/README.md).
