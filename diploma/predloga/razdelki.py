"""Dokonča izhod/diploma.docx po Vzorcu UL FKKT 2026 (kliče ga zgradi.sh po pandocu).

1. Odstrani Pandocov naslovni blok (Title/Author/Date).
2. Na začetek doda naslovnico, notranjo naslovnico in izjavo o avtorstvu iz metadata.yaml
   (besedilo izjave je prepisano iz Vzorca 2026, vključno z novo alinejo o UI).
3. Vsak del prvega odseka in vsako poglavje začne na novi LIHI strani (prelom razdelka »liha stran«),
   kot v Vzorcu 2026.
4. Prvi odsek (do Uvoda) nima številk strani; od Uvoda naprej so arabske številke spodaj na sredini, od 1.
5. Wordu naroči, naj ob odprtju posodobi polja (kazalo).
Uporaba: python3 predloga/razdelki.py izhod/diploma.docx metadata.yaml
"""
import copy
import sys
from pathlib import Path

import yaml
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm

docx_path, meta_path = Path(sys.argv[1]), Path(sys.argv[2])
meta = next(d for d in yaml.safe_load_all(meta_path.read_text(encoding="utf-8")) if d)
doc = Document(docx_path)
body = doc.element.body

# 1. Pandocov naslovni blok ven
for p in list(doc.paragraphs):
    if p.style is not None and p.style.name in ("Title", "Subtitle", "Author", "Date"):
        p._p.getparent().remove(p._p)

# 2. Naslovnici in izjava
spol = str(meta.get("spol") or "").strip().lower()
def g(m, z):
    """Moška / ženska oblika po polju spol v metadata.yaml; brez njega obe obliki kot v Vzorcu."""
    return m if spol == "m" else z if spol == "ž" else f"{m}/{z}"

ime = meta["author"]
naslov = meta["title"]
program = meta.get("program", "").strip()
program_vrstica = program if program.upper().startswith("UNIVERZITETNI") else f"UNIVERZITETNI ŠTUDIJSKI PROGRAM {program}"
kraj_leto = meta.get("kraj-leto", "Ljubljana")
mentor = meta.get("mentor", "")
mentor_rod = meta.get("mentor-rodilnik") or mentor
somentor = meta.get("somentor")

new = []  # (besedilo, slog, presledek_pred_pt)
new += [("UNIVERZA V LJUBLJANI", "Naslovnica", 0),
        ("FAKULTETA ZA KEMIJO IN KEMIJSKO TEHNOLOGIJO", "Naslovnica", 0),
        ("DIPLOMSKO DELO", "Naslovnica vrsta dela", 230),
        (ime, "Naslovnica", 30),
        (kraj_leto, "Naslovnica", 270)]
new += [("UNIVERZA V LJUBLJANI", "Naslovnica zacetek", 0),
        ("FAKULTETA ZA KEMIJO IN KEMIJSKO TEHNOLOGIJO", "Naslovnica", 0),
        (program_vrstica.upper(), "Naslovnica", 40),
        (naslov, "Naslovnica naslov dela", 90),
        ("DIPLOMSKO DELO", "Naslovnica vrsta dela", 40),
        (ime, "Naslovnica", 30),
        (f"MENTOR: {mentor}" if not mentor.lower().startswith("mentor") else mentor, "Naslovnica", 90)]
if somentor:
    new += [(f"SOMENTOR: {somentor}", "Naslovnica", 0)]
new += [(kraj_leto, "Naslovnica", 150 if not somentor else 130)]
alineje = [
    f"je diplomsko delo rezultat mojega raziskovalnega dela pod mentorstvom {mentor_rod}"
    + (f" in somentorstvom {meta.get('somentor-rodilnik') or somentor}" if somentor else "") + ";",
    f"sem {g('poskrbel', 'poskrbela')}, da so dela in mnenja drugih avtorjev, ki jih uporabljam v predloženem "
    "diplomskem delu, navedena oziroma citirana v skladu z navodili;",
    "se zavedam, da je plagiatorstvo, v katerem so tuje misli oziroma ideje predstavljene kot moje lastne, kaznivo "
    "po zakonu (Zakon o avtorski in sorodnih pravicah – uradno prečiščeno besedilo (ZASP-UPB3) (Ur. list RS, "
    "št. 16/2007));",
    f"sem {g('poskrbel', 'poskrbela')} za slovnično in oblikovno korektnost diplomskega dela;",
    "če bodo rezultati diplomskega dela po predhodnem soglasju mentorja"
    + (" in somentorja" if somentor else "") + " objavljeni, bo pri mojem imenu navedena afiliacija UL FKKT;",
    f"sem v primeru uporabe orodij umetne inteligence le-ta {g('uporabljal', 'uporabljala')} zgolj kot pomoč pri "
    "oblikovanju besedila ter zbiranju in obdelavi podatkov.",
]
new += [("IZJAVA O AVTORSTVU", "Izjava naslov", 0),
        ("diplomskega dela", "Izjava podnaslov", 0),
        (f"Spodaj {g('podpisani', 'podpisana')} {ime} sem {g('avtor', 'avtorica')} diplomskega dela z naslovom "
         f"{naslov}.", "Izjava besedilo", 36),
        ("S svojim podpisom zagotavljam, da:", "Izjava besedilo", 0)]
new += [("•\t" + a, "Izjava besedilo", 0) for a in alineje]
new += [(f"V Ljubljani, datum\t\tPodpis {g('avtorja', 'avtorice')}:", "Izjava besedilo", 48)]

first = body[0]
for text, style, before in new:
    p = doc.add_paragraph(text, style=style)
    if before:
        p.paragraph_format.space_before = Pt(before)
    if text.startswith("•\t"):
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
    first.addprevious(p._p)

# 3. + 4. Razdelki
ZACETKI = {"Naslovnica zacetek", "Izjava naslov", "Naslov povzetka", "Naslov brez stevilke", "Heading 1"}
final = body.find(qn("w:sectPr"))
paras = [p for p in doc.paragraphs]
starts = [p for i, p in enumerate(paras) if i > 0 and p.style is not None and p.style.name in ZACETKI]
first_main = next(p for p in starts if p.style.name == "Heading 1")
# del za vsakim razdelkom: ali je glavni (oštevilčen) in ali je prvi glavni
order = [paras[0]] + starts
is_main, seen = [], False
for p in order:
    seen = seen or (p is first_main)
    is_main.append(seen)


def make_sectpr(main, restart):
    sp = copy.deepcopy(final)
    for tag in ("w:headerReference", "w:footerReference", "w:titlePg", "w:pgNumType", "w:type"):
        for el in sp.findall(qn(tag)):
            sp.remove(el)
    if main:
        for el in final.findall(qn("w:footerReference")) + final.findall(qn("w:headerReference")):
            sp.insert(0, copy.deepcopy(el))
    t = OxmlElement("w:type"); t.set(qn("w:val"), "oddPage")
    pgsz = sp.find(qn("w:pgSz"))
    (pgsz.addprevious(t) if pgsz is not None else sp.append(t))
    if restart:
        n = OxmlElement("w:pgNumType"); n.set(qn("w:start"), "1")
        cols = sp.find(qn("w:cols"))
        (cols.addprevious(n) if cols is not None else sp.append(n))
    return sp


for idx in range(len(order) - 1):
    main = is_main[idx]
    restart = order[idx] is first_main
    brk = OxmlElement("w:p"); ppr = OxmlElement("w:pPr"); brk.append(ppr)
    ppr.append(make_sectpr(main, restart))
    order[idx + 1]._p.addprevious(brk)
# zadnji razdelek = končni sectPr
last = make_sectpr(is_main[-1], order[-1] is first_main)
final.getparent().replace(final, last)

# 5. Posodobi polja ob odprtju
settings = doc.settings.element
if settings.find(qn("w:updateFields")) is None:
    u = OxmlElement("w:updateFields"); u.set(qn("w:val"), "true")
    # vrstni red v settings.xml je predpisan (ECMA-376): updateFields stoji pred hdrShapeDefaults/footnotePr/...
    after = ("w:hdrShapeDefaults", "w:footnotePr", "w:endnotePr", "w:compat", "w:docVars", "w:rsids",
             "m:mathPr", "w:attachedSchema", "w:themeFontLang", "w:clrSchemeMapping",
             "w:doNotIncludeSubdocsInStats", "w:doNotAutoCompressPictures", "w:forceUpgrade", "w:captions",
             "w:readModeInkLockDown", "w:smartTagType", "w:schemaLibrary", "w:shapeDefaults",
             "w:doNotEmbedSmartTags", "w:decimalSymbol", "w:listSeparator")
    anchor = next((el for el in settings if any(el.tag == qn(t) for t in after)), None)
    (anchor.addprevious(u) if anchor is not None else settings.append(u))

doc.save(docx_path)
print(f"razdelki: {len(order)} (prvi odsek {is_main.count(False)}, drugi odsek {is_main.count(True)})")
