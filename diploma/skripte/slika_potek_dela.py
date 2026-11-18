"""Slika 3: shema poteka dela od surovih podatkov do vrednotenja napovedi.

Zagon iz korena repozitorija:
    python diploma/skripte/slika_potek_dela.py

Izhod:
    diploma/slike/potek_dela.png
"""
from __future__ import annotations

import os

KOREN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIKE = os.path.join(KOREN, 'diploma', 'slike')
CRNA = '#0b0b0b'
ROB = '#52514e'
POLNILO = '#f3f2ef'
MODALITETE = {'RNA': '#2a78d6', 'TCR': '#52514e', 'fluksi': '#eb6834'}


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'savefig.dpi': 300})

    s, v = 16.0, 8.2                       # cm
    fig = plt.figure(figsize=(s / 2.54, v / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, s); ax.set_ylim(0, v); ax.axis('off')

    def okvir(x, y, w, h, naslov, besedilo, barva=ROB):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.18',
                                    fc=POLNILO, ec=barva, lw=0.9))
        ax.text(x + w / 2, y + h - 0.17, naslov, ha='center', va='top', fontsize=8.4, fontweight='bold', color=CRNA)
        ax.text(x + w / 2, y + h - 0.6, besedilo, ha='center', va='top', fontsize=7.6, color=CRNA,
                linespacing=1.3)

    def puscica(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=9, lw=0.9,
                                     color=ROB, shrinkA=0, shrinkB=0))

    sredina = s / 2
    # 1 podatki in 2 kontrola kakovosti v isti vrsti
    y1, h1 = 6.3, 1.8
    okvir(0.2, y1, 7.3, h1, 'Podatki GSE200996 (3.1)',
          '28 bolnikov s HNSCC: štetja RNA,\nkontigi TCR in metapodatki celic T,\nkri in tumor, pred in po zdravljenju')
    puscica(7.5, y1 + h1 / 2, 8.5, y1 + h1 / 2)
    okvir(8.5, y1, 7.3, h1, 'Kontrola kakovosti (3.2)',
          '≥ 200 izraženih genov, celice iz metapodatkov,\n< 20 % mitohondrijskih štetij, veriga β s CDR3\n→ 69 930 celic')
    # 3 tri veje
    y3, h3, w3 = 3.9, 1.85, 5.0
    xs = [0.2, 5.5, 10.8]
    for x in xs:
        puscica(sredina, y1, x + w3 / 2, y3 + h3)
    okvir(xs[0], y3, w3, h3, 'Transkriptom (3.2)',
          'normalizirana štetja (enačba 3)\n↓\nPCA, 100 komponent', MODALITETE['RNA'])
    okvir(xs[1], y3, w3, h3, 'TCR (3.3)',
          'CDR3β → Atchleyjevi faktorji\n↓\nvložitev s 100 razsežnostmi',
          MODALITETE['TCR'])
    okvir(xs[2], y3, w3, h3, 'Metabolni fluksi (3.4)',
          'surova štetja 658 genov\n↓\nscFEA: 168 modulov (enačba 4)',
          MODALITETE['fluksi'])
    # 4 TRIM
    y4, h4 = 1.9, 1.45
    for x in xs:
        puscica(x + w3 / 2, y3, sredina + (x + w3 / 2 - sredina) * 0.55, y4 + h4)
    okvir(0.2, y4, 15.6, h4, 'Model TRIM (3.5–3.7)',
          'tri različice: RNA + TCR, fluksi + TCR, RNA + fluksi + TCR; izpusti enega bolnika (23 in 24); semena 0, 1, 2\n'
          'generiranje celic v štirih kvadrantih → psevdokloni (prag τ) → ocena $y_i$ (enačba 6)')
    puscica(sredina, y4, sredina, 1.5)
    # 5 vrednotenje
    okvir(0.2, 0.05, 15.6, 1.45, 'Vrednotenje (3.8 in 3.9)',
          'AUC z intervali bootstrap po klonih na celicah iz krvi pred zdravljenjem;\n'
          'primerjava z izhodiščnimi klasifikatorji (kNN, NN) in z velikostjo klona')

    os.makedirs(SLIKE, exist_ok=True)
    pot = os.path.join(SLIKE, 'potek_dela.png')
    fig.savefig(pot, facecolor='white', bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(pot)


if __name__ == '__main__':
    main()
