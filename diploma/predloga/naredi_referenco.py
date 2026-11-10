"""Iz privzete Pandoc reference.docx naredi predlogo z obliko po navodilih UL FKKT 2024:
A4, robovi 3 cm, Times New Roman 12, razmik 1,15, naslovi 16/14/12 (krepko / krepko / krepko ležeče),
številke strani spodaj na sredini. Uporaba: python3 predloga/naredi_referenco.py
"""
import subprocess
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

here = Path(__file__).parent
src = here / "_default.docx"
out = here / "fkkt_referenca.docx"
if not src.exists():
    subprocess.run(["pandoc", "-o", str(src), "--print-default-data-file", "reference.docx"], check=True)
doc = Document(src)

for s in doc.sections:
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    s.top_margin = s.bottom_margin = s.left_margin = s.right_margin = Cm(3)
    # številka strani spodaj na sredini
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
    if bold is not None: f.bold = bold
    if italic is not None: f.italic = italic

st = {x.name: x for x in doc.styles}
st.update({x.name.lower(): x for x in doc.styles})
for name in ("Normal", "Body Text", "First Paragraph", "Compact", "Bibliography", "Abstract",
             "Block Text", "Footnote Text", "Table", "Caption", "Image Caption", "Table Caption"):
    s = st.get(name) or st.get(name.lower())
    if s is None:
        continue
    font(s, 12)
    pf = s.paragraph_format
    pf.line_spacing = 1.15
    if name in ("Normal", "Body Text", "First Paragraph"):
        pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        pf.space_after = Pt(6)

for name, size, bold, italic in (("Heading 1", 16, True, False), ("Heading 2", 14, True, False),
                                 ("Heading 3", 12, True, True), ("Title", 16, True, False),
                                 ("Subtitle", 16, False, False), ("Author", 12, False, False),
                                 ("Date", 12, False, False)):
    s = st.get(name) or st.get(name.lower())
    if s is None:
        continue
    font(s, size, bold, italic)
    s.paragraph_format.space_before = Pt(18 if name == "Heading 1" else 12)
    s.paragraph_format.space_after = Pt(6)
    s.paragraph_format.line_spacing = 1.15
    if name == "Heading 1":
        s.paragraph_format.page_break_before = True
    if name in ("Title", "Subtitle", "Author", "Date"):
        s.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.save(out)
src.unlink()
print("shranjeno:", out)
