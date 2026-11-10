<!-- Cilj obsega: ~8–9 strani. Pregled literature, od splošnega k specifičnemu, konča z vrzeljo, ki jo zapolni naše delo.
     Vsaka trditev o biologiji/metodah mora imeti vir [@...]. Slike rišemo sami (»prirejeno po [n]«). -->

# Uvod

## Imunoterapija raka in limfociti T

<!-- - zaviralci imunskih kontrolnih točk (anti-PD-1), odziv le pri delu bolnikov
     - HNSCC: klinični kontekst, neoadjuvantna imunoterapija (kohorta GSE200996)
     - vloga CD8+ in CD4+ celic T v protitumorskem odzivu -->

## Receptor celic T in klonska ekspanzija

<!-- - V(D)J rekombinacija, raznolikost TCR, regija CDR3β kot glavni določevalec specifičnosti
     - klon = celice z istim TCR; klonska ekspanzija po aktivaciji
     - ekspanzija v krvi/tumorju po terapiji kot korelat odziva (literatura) -->

## Sekvenciranje posameznih celic: scRNA-seq in scTCR-seq

<!-- - 10x Genomics 5' + V(D)J, kaj dobimo (matrika štetij, kontigi)
     - kontrola kakovosti: prazne kapljice, delež mitohondrijskih genov, normalizacija
     - izzivi: redkost (dropout), paketni učinki, multipleksiranje -->

## Globoki generativni modeli za podatke posameznih celic

<!-- - avtokodirnik in VAE (ELBO, latentni prostor), scVI kot zgled
     - multimodalni VAE (product/mixture of experts, povprečenje latentnih porazdelitev)
     - TRIM: arhitektura, pogojevanje na tkivo in fazo zdravljenja, generiranje celic za napoved ekspanzije, psevdokloni -->

## Celični metabolizem limfocitov T

<!-- - metabolno preprogramiranje ob aktivaciji (glikoliza, OXPHOS, glutaminoliza), naivne vs efektorske vs spominske celice
     - metabolizem v tumorskem mikrookolju (hipoksija, pomanjkanje glukoze), izčrpanost T-celic
     - zakaj bi metabolno stanje lahko napovedovalo proliferacijo/ekspanzijo -->

## Ocenjevanje metabolnih fluksov iz transkriptomskih podatkov

<!-- - genomski metabolni modeli (Human-GEM), stehiometrična matrika, FBA, pFBA
     - kontekstno specifični modeli (iMAT, INIT/ftINIT, RAVEN)
     - metode za posamezne celice: Compass, scFEA, METAFlux — prednosti in omejitve
     - pravila GPR in ocene aktivnosti reakcij -->

## Vrzel v znanju

<!-- 1 odstavek: TRIM uporablja le RNA + TCR; metabolna informacija ni eksplicitno vključena;
     ni znano, ali fluksi kot dodatna modaliteta izboljšajo napoved klonske ekspanzije. -->
