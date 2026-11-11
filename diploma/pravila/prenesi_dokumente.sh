#!/usr/bin/env bash
# Prenese uradne dokumente (UL FKKT + UL) za diplomsko delo v mapo ./originali/
# Zaženi lokalno:  bash diploma/pravila/prenesi_dokumente.sh
set -u
cd "$(dirname "$0")"
mkdir -p originali
dl() { echo "-> $1"; curl -fsSL --retry 3 -o "originali/$1" "$2" && return 0; echo "   NI USPELO: $2"; return 1; }

# --- UL FKKT: pravilniki, navodila in vzorec (veljajo od 1. 10. 2026; povezave preverjene 11. 11. 2026) ---
B=https://fkkt.uni-lj.si/fileadmin/user_upload
dl 01_FKKT_Pravilnik_o_diplomskem_delu_preciscen_2026.pdf "$B/habilitacije/PRAVILNIK_O_DIPLOMSKEM_DELU_PRVOSTOPENJSKIH_STUDIJSKIH_PROGRAMOV_FAKULTETE_ZA_KEMIJO_IN_KEMIJSKO_TEHNOLOGIJO_UNIVERZE_V_LJUBLJANI.pdf"
dl 01a_FKKT_Pravilnik_o_spremembah_2026.pdf "$B/o_fakulteti/Pravilnik_o_spremembah_Pravilnika_o_diplomskem_delu_prvostopenjskih_studijskih_programov_Fakultete_za_kemijo_in_kemijsko_tehnologijo_Univerze_v_Ljubl.pdf"
dl 02_FKKT_Navodila_za_izdelavo_dipl_mag_dela_2026.docx "$B/habilitacije/Navodila_za_izdelavo_dipl_mag_dela_2026_-_cistopis.docx"
dl 03_FKKT_Rokovnik_dejavnosti_za_diplomsko_delo.pdf "$B/Studentski_referat/Pravilniki_studij/ROKOVNIK_DEJAVNOSTI_ZA_DIPLOMSKO_DELO.pdf"
dl 10_FKKT_Vzorec_diploma_magisterij_2026.docx "$B/habilitacije/Vzorec_diploma_magisterij_2026_-_cistopis.docx"
dl 11_FKKT_Template_diploma_master_2026_EN.docx "$B/Studentski_referat/Pravilniki_studij/Template_diploma_master_2026.docx"

# --- Primer: F. Petrovič, diplomsko delo UL FKKT 2025 (Biokemija, mentor M. Moškon) ---
dl 30_Primer_Petrovic_2025_ftINIT.pdf "https://repozitorij.uni-lj.si/Dokument.php?id=215360&lang=slv"

# --- Wiki FKKT (Biokemija): obrazec za seznam nebesedilnih sestavin (strežnik je včasih nedosegljiv) ---
dl 12_FKKT_Biokemija_Seznam_nebesedilnih_sestavin_2022c.xlsx "http://web.fkkt.uni-lj.si/biokemija/Seznam%20nebesedilnih%20sestavin%202022c.xlsx"

# --- UL: umetna inteligenca ---
dl 20_UL_in_UI_izjava_in_smernice_2023.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/UL_in_UI_izjava_in_smernice_za_objavo_1_.pdf"
dl 21_UL_Smernice_za_uporabo_UI_2026.pdf "https://www.teof.uni-lj.si/uploads/UL26_smernice-UI(1).pdf"
dl 22_UL_Odgovorno_in_eticno_raziskovanje.pdf "https://www.uni-lj.si/assets/Sluzba-za-raziskovalno-dejavnost/Etika-in-integriteta/Odgovorno-in-eticno-raziskovanje.pdf"

# --- Arhiv (starejše verzije, samo za primerjavo; NE uporabljaj) ---
dl 90_ARHIV_FKKT_Navodila_izdelava_diploma_magisterij_koncna.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/Studij/pravilniki/Navodila_izdelava_diploma_magisterij_koncna.pdf"
dl 91_ARHIV_FKKT_Navodila_za_izdelavo_dipl_mag_dela_2018.pdf "https://fkkt.uni-lj.si/fileadmin/user_upload/arhivske_datoteke/1-O_fakulteti/3-Pravilniki__akti/2_stopnja_Pravilnik_Priloge/Navodila_za_izdelavo_dipl_mag_dela_2018.pdf"

echo "Končano. Datoteke so v $(pwd)/originali"
