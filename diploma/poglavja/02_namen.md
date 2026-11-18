<!-- Navodila Biokemija v1.2 + kontrolni seznam:
     - naslov po številu hipotez (tri hipoteze → »Namen dela in hipoteze«);
     - pisano v PRETEKLIKU (kaj smo želeli ugotoviti pred začetkom dela);
     - namen kot opis problema, ne povzetek dela; lahko predvidene metode, NE rezultati;
     - vsaj ena hipoteza kot trditev, katere resničnosti pred delom ni bilo mogoče vedeti. -->

# Namen dela in hipoteze

Delitev celic T je metabolno zahtevna (razdelek 1.6), zato smo predvidevali, da bi metabolno stanje celice T pred zdravljenjem lahko vsebovalo informacijo o tem, ali se bo njen klon po sprostitvi zavore PD-1 namnožil. Namen dela je bil ugotoviti, ali metabolno stanje celic T, ocenjeno iz podatkov scRNA-seq, takšno informacijo vsebuje in ali jo model TRIM lahko izkoristi. Zanimalo nas je predvsem, ali opis celice z metabolnimi fluksi napovedi doda novo informacijo ali pa le drugače zapiše informacijo, ki jo transkriptom že vsebuje.

Delo smo načrtovali v treh korakih. Najprej smo želeli model TRIM ponoviti na javno dostopnih podatkih bolnikov s HNSCC, zdravljenih z neoadjuvantno imunoterapijo [@luoma2022], in pri tem podatke predobdelati sami. Nato smo nameravali z metodo scFEA [@alghamdi2021] za vsako celico oceniti metabolne flukse skozi metabolne module in jih v model TRIM vključiti kot tretjo modaliteto: v eni različici namesto transkriptoma, v drugi poleg njega. Na koncu smo želeli vse različice ovrednotiti z isto shemo »izpusti enega bolnika«, jih primerjati z izhodiščnimi klasifikatorji in med seboj ter vsako učenje ponoviti z več semeni generatorja naključnih števil.

Postavili smo tri hipoteze:

- **H1:** Model TRIM, ponovljen na javnih podatkih z lastno predobdelavo, napove klonsko ekspanzijo v krvi izpuščenega bolnika z višjo AUC kot izhodiščni klasifikatorji na istih celicah.
- **H2:** Metabolni fluksi, ocenjeni s scFEA, skupaj z receptorjem celic T vsebujejo informacijo o klonski ekspanziji: model TRIM s fluksi in TCR napove ekspanzijo z višjo AUC kot izhodiščni klasifikatorji.
- **H3:** Dodajanje metabolnih fluksov k transkriptomu in receptorju celic T izboljša napoved klonske ekspanzije v primerjavi z modelom TRIM, ki uporablja samo transkriptom in TCR.
