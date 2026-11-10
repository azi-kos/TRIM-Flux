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
├── predloga/              Word referenca z obliko FKKT (A4, 3 cm, TNR 12, 1,15, naslovi 16/14/12)
├── pravila/               pravila FKKT/UL, citatni slog fkkt-ul.csl, orodja za literaturo
└── zgradi.sh              Markdown -> izhod/diploma.docx
```

## Gradnja

```bash
bash diploma/zgradi.sh                 # potrebuje pandoc (>= 3) in python-docx
```

Kazalo se v Wordu osveži z desnim klikom → »Posodobi polje«. Naslovnico, notranjo naslovnico in Izjavo o avtorstvu na koncu prenesemo iz uradnega vzorca FKKT (`pravila/prenesi_dokumente.sh` ga prenese). Sledi izvoz v PDF/A.

## Pisanje

- Citat: `... klonsko ekspanzijo [@kljuc2023].` V Wordu postane `[1]`, pred piko, oštevilčeno po vrstnem redu.
- Komentarji `<!-- ... -->` so delovne opombe in se ne izpišejo.
- Obseg: okoli 30 strani. Okvirno: Uvod 8–9, Namen 1–2, Metode 7–8, Rezultati in razprava 8–9, Zaključek 1–2.
- Številke v rezultatih vpisujemo iz shranjenih logov v `results/logs/`, ne na pamet.
- Pravila oblike, citiranja in uporabe UI so povzeta v [`pravila/README.md`](pravila/README.md).
