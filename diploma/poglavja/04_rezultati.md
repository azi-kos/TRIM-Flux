<!-- ~8–9 strani. Rezultati s tabelami in slikami (lastne, iz results/ in analize/), razprava primerja s hipotezami in literaturo.
     Vse številke vpisujemo iz shranjenih JSON logov (results/logs/), ne na pamet. -->

# Rezultati in razprava

## Kontrola kakovosti in opis podatkov

<!-- - SLIKA: porazdelitev št. genov na celico (bimodalnost, prag 200), delež mito
     - TABELA: celice po bolnikih/kvadrantih po QC -->

## Ponovitev modela TRIM (H1)

<!-- - vpliv normalizacije: log1p → normalizacija avtorjev (AUC P24 0,77 → 0,87); izhodišče knn_both 0,66
     - SLIKA: krivulja ROC, UMAP prekrivanje generiranih in resničnih celic
     - primerjava z objavljenimi rezultati TRIM -->

## Metabolni fluksi kot samostojna modaliteta (H2)

<!-- - flux_tcr: AUC (Monte Carlo, povprečje ± std po semenih, P23/P24 ...)
     - interpretacija: nad izhodiščem, pod RNA + TCR -->

## Združevanje RNA, TCR in fluksov (H3)

<!-- - rna_flux_tcr vs rna_tcr; enkratni 0,954 je bil šum (pokaži razpršenost med semeni)
     - statistična primerjava, interval zaupanja; ablacija (var. 3), če bo izvedena -->

## Razprava

<!-- - zakaj fluksi ne prispevajo (veliko) nove informacije: izpeljani iz iste RNA, kompresija
     - omejitve: malo bolnikov, nedeterminizem učenja, izbira metode za fluks, kalibracija praga
     - primerjava z literaturo o metabolizmu T-celic -->
