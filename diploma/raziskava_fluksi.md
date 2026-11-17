# Katera metoda za flukse celic T iz scRNA-seq (17. 11. 2026)

Vprašanje: katera metoda najbolje oceni metabolne flukse celic T iz podatkov RNA, če potrebujemo vrednosti za vsako od ~70 000 celic (vhod v TRIM)?

## Kratek odgovor

Nobena metoda ni preverjena na ravni fluksov posamezne celice. Vse preverbe so na ravni populacij celic. Ocene po celicah se med metodami skoraj ne ujemajo: korelacija med scFEA in Compass po celicah je od −0,26 do 0,04 [1]. Izbira je zato odvisna od tega, kaj metoda zmore pri naših podatkih.

- **Najbolje preverjena pri celicah T: Compass** [2]. Preverjena je na celicah Th17 (LC-MS, sledenje ¹³C, Seahorse, farmakološke in genske motnje, EAE in vivo). Ni primerna za 70 000 posameznih celic: dokumentacija priporoča psevdobulk, Turbo-Compass na velikem naboru še vedno teče ~36 h s 50 procesorji [3]. Izhod so ocene reakcij (kazni), ne fluksi.
- **Najboljša v primerjavi 2026 za modele posameznih celic: ftINIT** [4] (z normalizacijo in imputacijo, prag 25. percentil). Na drugem mestu je GIMME, ki je tudi najhitrejši. ftINIT je počasnejši in ima časovne prekinitve, zato modela za vsako od 70 000 celic ni realno zgraditi. Gustafsson in sod. [5] pokažejo, da so za stabilen model iz kapljičnih podatkov potrebne tisoče celic (preverjeno tudi na celicah T). ftINIT potrebuje Gurobi.
- **METAFlux** [6] flukse za scRNA-seq računa samo na ravni celičnih tipov ali gruč: »Here, we do not compute cell-wise flux« [7].
- **scCellFie** [8] hitro izračuna aktivnost 218 metabolnih nalog za vsako celico (30 milijonov celic), a to niso fluksi. Ni eksperimentalne preverbe, je preprint.
- **scFEA** [9] (dosedanja metoda) da vrednosti po celicah, a je stohastičen, deluje na modulih in se po celicah ne ujema s Compass [1]. Študija CD8 T [10] ga je uporabila brez lastne preverbe.

## Priporočilo za diplomo

**ftINIT csGEM za celice T in nato fluksi po celicah s kvadratnim programom v tem modelu.**

1. En csGEM iz Human-GEM z ftINIT na združenih celicah T vseh bolnikov razen 23 in 24. Bazen ima več deset tisoč celic, kar zadosti [5], in izpuščena bolnika ne vplivata na mrežo.
2. Za vsako celico meje reakcij iz izražanja in pravil GPR v csGEM, nato strogo konveksni kvadratni program (kot METAFlux v načinu za bulk [6]). Rezultat je enoličen in determinističen. Redkost podatkov ublažimo z glajenjem med sosednjimi celicami istega vzorca (bolnik, tkivo, čas) ali z imputacijo, kot priporoča [4].
3. Reakcije skrčimo (podsistemi ali najbolj variabilne reakcije) za vhod v TRIM.

Rezervna pot: isti kvadratni program na celotnem Human-GEM (brez ftINIT), če Gurobi ali ftINIT ne steče. Kodo za korak 2 zato najprej preizkusimo na Human-GEM.

Zakaj: izhod so pravi fluksi (ustreza prijavljenemu naslovu »metabolnih fluksov«), mreža je prilagojena celicam T, rezultat je ponovljiv brez semen in kosov, ftINIT pa je metoda, ki jo mentor pozna.

## Viri

1. Amiss JMB, Jabbari H, Lum J. A method to identify high consensus predictions of single-cell metabolic flux. bioRxiv 2024. https://doi.org/10.1101/2024.01.15.572211
2. Wagner A in sod. Metabolic modeling of single Th17 cells reveals regulators of autoimmunity. Cell 2021. https://doi.org/10.1016/j.cell.2021.05.045 (v bib: wagner2021)
3. Dokumentacija Compass, Turbo-Compass in nastavitve: https://compass-wagnerlab.readthedocs.io/en/latest/turbo_compass.html
4. Yang J, Deng Y, Luo J, Wang Y, Li F, Chen Y. A benchmarking framework for single-cell genome-scale metabolic model construction. bioRxiv 2026. https://doi.org/10.64898/2026.08.02.742346
5. Gustafsson J in sod. Generation and analysis of context-specific genome-scale metabolic models derived from single-cell RNA-Seq data. PNAS 2023;120(6):e2217868120. https://doi.org/10.1073/pnas.2217868120 (v bib: gustafsson2023)
6. Huang Y in sod. Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux. Nat Commun 2023;14:4883. https://doi.org/10.1038/s41467-023-40457-w (v bib: huang2023)
7. Vodič METAFlux, razdelek za scRNA-seq: https://github.com/KChen-lab/METAFlux (Tutorials/pipeline.Rmd)
8. Armingol E in sod. Atlas-scale metabolic activities inferred from single-cell and spatial transcriptomics (scCellFie). bioRxiv 2025. https://doi.org/10.1101/2025.05.09.653038
9. Alghamdi N in sod. A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data. Genome Res 2021;31(10):1867–1884. https://doi.org/10.1101/gr.271205.120 (v bib: alghamdi2021)
10. In silico single-cell metabolism analysis unravels a new transition stage of CD8 T cells 4 days post-infection. bioRxiv 2023 / ImmunoInformatics 2024. https://www.biorxiv.org/content/10.1101/2023.09.22.558248v2
