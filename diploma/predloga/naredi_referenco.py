"""Iz privzete Pandoc reference.docx naredi predlogo z obliko po Navodilih in Vzorcu UL FKKT 2026.

Vrednosti so prepisane iz Vzorec_diploma_magisterij_2026_-_cistopis.docx (styles.xml, prebrano 11. 11. 2026)
in Navodil 2026:
- A4, robovi 3 cm (1701 twip), glava/noga 2 cm;
- besedilo: Times New Roman 12, obojestransko, 12 pt za odstavkom, razmik 1,15 (Navodila; Vzorec ima 1,3);
- naslov 1. ravni: 16 krepko, levo, 30 pt za njim; 2. raven: 14 krepko, 24 pt pred; 3. raven: 12 krepko ležeče,
  18 pt pred, 6 pt za;
- napisi slik/tabel: 12, sredinsko, enojni razmik, 6 pt pred, 24 pt za;
- neoštevilčeni naslovi (Kazalo, Seznam uporabljenih kratic in simbolov): 16 krepko kot Heading 5 v Vzorcu;
- številka strani spodaj na sredini (nogo ima samo drugi odsek, to uredi razdelki.py).
Uporaba: python3 predloga/naredi_referenco.py
"""
import subprocess
from pathlib import Path
from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm, RGBColor

here = Path(__file__).parent
src = here / "_default.docx"
out = here / "fkkt_referenca.docx"
if not src.exists():
    subprocess.run(["pandoc", "-o", str(src), "--print-default-data-file", "reference.docx"], check=True)
doc = Document(src)

for s in doc.sections:
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    s.top_margin = s.bottom_margin = s.left_margin = s.right_margin = Cm(3)
    s.header_distance = s.footer_distance = Cm(2)
    p = s.footer.paragraphs[0] if s.footer.paragraphs else s.footer.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    for tag, text in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), tag)
        else:
            el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = text
        run._r.append(el)


def font(style, size, bold=None, italic=None):
    f = style.font
    f.name = "Times New Roman"
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts"); rpr.append(rfonts)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rfonts.set(qn(a), "Times New Roman")
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rfonts.get(qn(a)) is not None:
            del rfonts.attrib[qn(a)]
    f.size = Pt(size)
    f.color.rgb = RGBColor(0, 0, 0)
    if bold is not None:
        f.bold = bold
    if italic is not None:
        f.italic = italic


def para(style, before=None, after=None, line=None, align=None, keep=None):
    pf = style.paragraph_format
    if before is not None:
        pf.space_before = Pt(before)
    if after is not None:
        pf.space_after = Pt(after)
    if line is not None:
        pf.line_spacing = line
    if align is not None:
        pf.alignment = align
    if keep is not None:
        pf.keep_with_next = keep


st = {x.name.lower(): x for x in doc.styles}


def get_or_add(name, base="Normal"):
    s = st.get(name.lower())
    if s is None:
        s = doc.styles.add_style(name, WD_STYLE_TYPE.PARAGRAPH)
        s.base_style = st["normal"]
        st[name.lower()] = s
    return s


J, C, L = WD_ALIGN_PARAGRAPH.JUSTIFY, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.LEFT

for name in ("Normal", "Body Text", "First Paragraph", "Compact", "Bibliography", "Abstract",
             "Block Text", "Footnote Text", "Table"):
    s = st.get(name.lower())
    if s is None:
        continue
    font(s, 12)
    para(s, line=1.15)
    if name in ("Normal", "Body Text", "First Paragraph"):
        para(s, before=0, after=12, align=J)
    if name == "Bibliography":
        para(s, before=0, after=6, align=L)

for name in ("Caption", "Image Caption", "Table Caption"):
    s = st.get(name.lower())
    if s is None:
        continue
    font(s, 12, bold=False, italic=False)
    para(s, before=6, after=24, line=1.0, align=C)

for name, size, bold, italic, before, after in (
        ("Heading 1", 16, True, False, 0, 30),
        ("Heading 2", 14, True, False, 24, 0),
        ("Heading 3", 12, True, True, 18, 6)):
    s = st.get(name.lower())
    font(s, size, bold, italic)
    para(s, before=before, after=after, line=1.15, align=L, keep=True)
    s.paragraph_format.page_break_before = False  # nove strani (lihe) doda razdelki.py

# Neoštevilčeni naslovi prvega odseka (kot Heading 5 v Vzorcu): niso v kazalu.
s = get_or_add("Naslov brez stevilke")
font(s, 16, True, False)
para(s, before=0, after=30, line=1.15, align=L, keep=True)

# Naslov dela na strani s povzetkom (Navodila 2026, točka 2.7: velikost 12, krepko).
s = get_or_add("Naslov povzetka")
font(s, 12, True, False)
para(s, before=0, after=12, line=1.15, align=L, keep=True)

# Naslovnici in izjava o avtorstvu (Vzorec 2026: 16 pt sredinsko, »DIPLOMSKO DELO« 20 pt, izjava 18/16).
for name, size, bold in (("Naslovnica", 16, False), ("Naslovnica zacetek", 16, False),
                         ("Naslovnica naslov dela", 16, True), ("Naslovnica vrsta dela", 20, False),
                         ("Izjava naslov", 18, True), ("Izjava podnaslov", 16, False)):
    s = get_or_add(name)
    font(s, size, bold, False)
    para(s, before=0, after=0, line=1.15, align=C)
s = get_or_add("Izjava besedilo")
font(s, 12, False, False)
para(s, before=0, after=12, line=1.15, align=L)

doc.save(out)
src.unlink()
print("shranjeno:", out)
