Table: Tabela 4: AUC napovedi klonske ekspanzije modela TRIM pri izpuščenih bolnikih 23 in 24. V oklepaju je 95-% interval zaupanja po DeLongu; »–« pomeni, da zagona nismo izvedli. Povprečje je izračunano po razpoložljivih semenih.

| Različica | Seme | Bolnik 23 | Bolnik 24 |
|----------------------|-------------|--------------------------|--------------------------|
| RNA + TCR | 0 | 0,906 (0,887–0,925) | 0,916 (0,856–0,976) |
|  | 1 | 0,782 (0,754–0,809) | 0,773 (0,671–0,876) |
|  | 2 | 0,923 (0,905–0,941) | – |
|  | povprečje | 0,870 | 0,845 |
| fluksi + TCR | 0 | 0,745 (0,715–0,774) | 0,813 (0,720–0,906) |
|  | 1 | 0,476 (0,442–0,510) | 0,659 (0,567–0,750) |
|  | 2 | 0,780 (0,752–0,808) | 0,851 (0,779–0,922) |
|  | povprečje | 0,667 | 0,774 |
| RNA + fluksi + TCR | 0 | 0,742 (0,713–0,771) | 0,848 (0,757–0,939) |
|  | 1 | 0,908 (0,889–0,928) | 0,806 (0,706–0,907) |
|  | 2 | 0,909 (0,890–0,928) | 0,785 (0,681–0,889) |
|  | povprečje | 0,853 | 0,813 |

Table: Tabela 5: AUC izhodiščnih klasifikatorjev. Pri vhodih, ki smo jih uporabili v več zagonih, je naveden razpon vrednosti.

| Klasifikator | Vhod | Bolnik 23 | Bolnik 24 |
|---|---|---|---|
| kNN | RNA | 0,485 | 0,652 |
| kNN | TCR | 0,518 | 0,419 |
| kNN | RNA in TCR | 0,488 | 0,659 |
| kNN | fluksi | 0,506–0,515 | 0,461–0,486 |
| kNN | fluksi in TCR | 0,464 | 0,206 |
| NN | RNA | 0,472 | 0,348 |
| NN | TCR | 0,476–0,481 | 0,457–0,477 |
| NN | RNA in TCR | 0,397–0,444 | 0,153–0,191 |
| NN | fluksi | 0,541–0,548 | 0,437–0,496 |
| NN | fluksi in TCR | 0,628–0,733 | 0,206–0,503 |

Table: Tabela 6: Razlika v AUC (ΔAUC) različic RNA + fluksi + TCR in fluksi + TCR do različice RNA + TCR na istih celicah pri istem bolniku in semenu, s parnim DeLongovim testom. Pozitivna razlika pomeni višjo AUC različice s fluksi; pri bolniku 24 s semenom 2 zagona RNA + TCR nismo izvedli.

| Bolnik | Seme | RNA + fluksi + TCR: ΔAUC | *p* | fluksi + TCR: ΔAUC | *p* |
|----------|--------|------------------------|------------|------------------------|------------|
| 23 | 0 | −0,164 | < 0,001 | −0,161 | < 0,001 |
| 23 | 1 | 0,127 | < 0,001 | −0,306 | < 0,001 |
| 23 | 2 | −0,014 | 0,17 | −0,143 | < 0,001 |
| 24 | 0 | −0,068 | 0,18 | −0,104 | 0,014 |
| 24 | 1 | 0,033 | 0,45 | −0,114 | 0,020 |
