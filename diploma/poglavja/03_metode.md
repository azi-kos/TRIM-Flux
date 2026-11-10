<!-- ~7–8 strani. Dovolj natančno za ponovitev; logični (ne kronološki) vrstni red.
     Vsaka metoda iz literature z virom. Navedi verzije programske opreme in strojno opremo. -->

# Materiali in metode

## Podatki

<!-- - GEO GSE200996: HNSCC, neoadjuvantna imunoterapija; kri in tumor, pred in po zdravljenju
     - vsebina: .h5 matrike štetij, filtered_contig (TCR), metapodatki po celicah
     - število bolnikov (28), celic T po QC, porazdelitev po kvadrantih (tkivo × faza)
     - TABELA: št. celic po bolniku in kvadrantu (iz notebooks_flux/00_diagnostika_pacientov.ipynb) -->

## Predobdelava

<!-- - nalaganje (scanpy), ujemanje črtnih kod s ključem (barcode, bolnik) zaradi multipleksiranja
     - kontrola kakovosti: min_genes = 200 (bimodalna porazdelitev), delež mito < 20 %
     - normalizacija po avtorjih TRIM (log, odštej min, deli z max na celico) – enačba
     - TCR: izbira kontiga z največ UMI, CDR3β, V in J geni; štetje klonov po kvadrantih
     - PCA (100 komponent) na RNA
     - znana omejitev: ~0,3 % celic z napačno oznako tkiva/faze zaradi kolizij ključev -->

## Vložitev zaporedij TCR

<!-- - konvolucijski avtokodirnik na CDR3β (learn_tcr_embedding.py), 100-dimenzionalna vložitev, parametri učenja -->

## Model TRIM

<!-- - SLIKA (lastna): arhitektura – kodirniki po modalitetah, povprečenje latentnih porazdelitev, pogojni dekodirniki
     - funkcija izgube (rekonstrukcija RNA, TCR, KL, uteži λ), hiperparametri (20 000 korakov ...)
     - psevdokloni: L1 razdalje med generiranimi TCR, prag kalibriran na krvi pred zdravljenjem
     - prilagoditve za omejen RAM (razdalje samo na izpuščenem bolniku) in razlika od originala pri kalibraciji praga -->

## Izračun metabolnih fluksov

<!-- !!! Odvisno od končne izbire metode (scFEA že izvedeno; načrtovan prehod na pristop na osnovi GEM/GPR).
     - scFEA: moduli (~168), vhod surova štetja, učenje na GPU, pospešitve (vektorizacija)
     - [načrtovano] Human-GEM / csGEM (ftINIT), ocene aktivnosti reakcij prek GPR ali QP po celici
     - standardizacija z (povprečje 0, std 1 po modulu) – zakaj (skala izgube) -->

## Razširitev modela TRIM s fluksi

<!-- - variante: rna_tcr (izvirnik), flux_tcr (var. 1), rna_flux_tcr (var. 2), var. 3 (RNA brez genov fluksov + fluksi + TCR)
     - kodirnik/dekodirnik za fluks, uteži izgube -->

## Vrednotenje napovedi klonske ekspanzije

<!-- - shema »izpusti enega bolnika«; izbira bolnikov (vsi 4 kvadranti, dovolj ekspandiranih klonov)
     - resnica: klon ekspandiran, če blood-post > blood-pre; napoved: generirani psevdokloni
     - metrike: ROC AUC, priklic/natančnost pri Youdenovem pragu; Hanley–McNeil CI; več semen (Monte Carlo)
     - izhodiščni klasifikatorji (kNN, NN na RNA/TCR/obeh) -->

## Programska in strojna oprema

<!-- Python, PyTorch, scanpy, scikit-learn (verzije iz environment.yml); Google Colab, NVIDIA A100 80 GB.
     Koda: https://github.com/azi-kos/TRIM-Flux (veja flux-extension). -->
