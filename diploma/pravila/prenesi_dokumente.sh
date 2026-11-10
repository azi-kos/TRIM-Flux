#!/usr/bin/env bash
# Prenese uradne dokumente (UL FKKT + UL) za diplomsko delo v mapo ./originali/
# Zaženi lokalno:  bash diploma/pravila/prenesi_dokumente.sh
set -u
cd "$(dirname "$0")"
mkdir -p originali
dl() { echo "-> $1"; curl -fsSL --retry 3 -o "originali/$1" "$2" || echo "   NI USPELO: $2"; }

# --- UL FKKT: pravilniki in navodila ---
dl 01_FKKT_Pravilnik_o_diplomskem_delu_cistopis_okt2023.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/cistopis_PRAVILNIK_O_DIPLOMSKEM_DELU_okt_2023.pdf"
dl 02_FKKT_Navodila_za_izdelavo_dipl_mag_dela_2024.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki_in_obrazci/Navodila_za_izdelavo_dipl_mag_dela_2024.pdf"
dl 03_FKKT_Rokovnik_dejavnosti_za_diplomsko_delo.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/ROKOVNIK_DEJAVNOSTI_ZA_DIPLOMSKO_DELO.pdf"
dl 04_FKKT_Pravilnik_1_stopnja.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/Pravilnik_1_stopnja.pdf"

# --- UL FKKT: Word predlogi (uradni vzorec) ---
dl 10_FKKT_Vzorec_diploma_magisterij_2024.docx "https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Vzorec_diploma_magisterij_2024.docx"
dl 11_FKKT_Template_diploma_master_2024_EN.docx "https://fkkt.uni-lj.si/fileadmin/user_upload/Studentski_referat/Pravilniki_studij/Template_diploma_master_2024.docx"

# --- UL: umetna inteligenca ---
dl 20_UL_in_UI_izjava_in_smernice_2023.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/UL_in_UI_izjava_in_smernice_za_objavo_1_.pdf"
dl 21_UL_Smernice_za_uporabo_UI_2026.pdf "https://www.teof.uni-lj.si/uploads/UL26_smernice-UI(1).pdf"
dl 22_UL_Odgovorno_in_eticno_raziskovanje.pdf "https://www.uni-lj.si/assets/Sluzba-za-raziskovalno-dejavnost/Etika-in-integriteta/Odgovorno-in-eticno-raziskovanje.pdf"

# --- Arhiv (starejše verzije, samo za primerjavo; NE uporabljaj) ---
dl 90_ARHIV_FKKT_Navodila_izdelava_diploma_magisterij_koncna.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/Navodila_izdelava_diploma_magisterij_koncna.pdf"
dl 91_ARHIV_FKKT_Navodila_za_izdelavo_dipl_mag_dela_2018.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/arhivske_datoteke/1-O_fakulteti/3-Pravilniki__akti/2_stopnja_Pravilnik_Priloge/Navodila_za_izdelavo_dipl_mag_dela_2018.pdf"

echo "Končano. Datoteke so v $(pwd)/originali"
