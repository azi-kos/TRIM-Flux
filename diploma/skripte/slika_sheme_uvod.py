"""Sliki 1 in 2 v Uvodu: shema modela TRIM s postopkom napovedi in shema metode scFEA.

Obe sta lastni poenostavljeni shemi, prirejeni po He in sod. (2026) oziroma Alghamdi in sod. (2021).

Zagon iz korena repozitorija:
    python diploma/skripte/slika_sheme_uvod.py

Izhod:
    diploma/slike/shema_trim.png
    diploma/slike/shema_scfea.png
"""
from __future__ import annotations

import os

KOREN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SLIKE = os.path.join(KOREN, 'diploma', 'slike')
CRNA = '#0b0b0b'
ROB = '#52514e'
SIVA = '#a3a29b'
POLNILO = '#f3f2ef'
RNA = '#2a78d6'
TCR = '#52514e'
FLUKS = '#eb6834'


def pripravi_plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'savefig.dpi': 300, 'mathtext.fontset': 'dejavusans'})
    return plt


def platno(plt, s, v):
    fig = plt.figure(figsize=(s / 2.54, v / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, s)
    ax.set_ylim(0, v)
    ax.axis('off')
    return fig, ax


def okvir(ax, x, y, w, h, besedilo, barva=ROB, fs=7.5, krepko=False, polnilo=POLNILO, lw=0.9, ls='-'):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.15',
                                fc=polnilo, ec=barva, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, besedilo, ha='center', va='center', fontsize=fs, color=CRNA,
            fontweight='bold' if krepko else 'normal', linespacing=1.3)


def trapez(ax, x, y, w, h, besedilo, barva, ozek_desno=True):
    from matplotlib.patches import Polygon
    z = 0.28 * h
    if ozek_desno:
        tocke = [(x, y), (x + w, y + z), (x + w, y + h - z), (x, y + h)]
    else:
        tocke = [(x, y + z), (x + w, y), (x + w, y + h), (x, y + h - z)]
    ax.add_patch(Polygon(tocke, closed=True, fc=POLNILO, ec=barva, lw=0.9))
    ax.text(x + w / 2, y + h / 2, besedilo, ha='center', va='center', fontsize=7.5, color=CRNA, linespacing=1.3)


def puscica(ax, x0, y0, x1, y1, barva=ROB, ls='-', lw=0.9, rad=0.0):
    from matplotlib.patches import FancyArrowPatch
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, lw=lw, color=barva,
                                 linestyle=ls, shrinkA=0, shrinkB=0,
                                 connectionstyle=f'arc3,rad={rad}'))


def crka(ax, x, y, c, besedilo):
    ax.text(x, y, c, fontsize=10, fontweight='bold', va='center', ha='left', color=CRNA)
    ax.text(x + 0.45, y, besedilo, fontsize=8.5, va='center', ha='left', color=CRNA)


def shema_trim(plt):
    import numpy as np
    s, v = 16.0, 10.6
    fig, ax = platno(plt, s, v)

    # A: zgradba modela
    crka(ax, 0.1, 10.25, 'A', 'zgradba modela in učenje')
    yr, yt, h = 7.75, 5.55, 1.15
    okvir(ax, 0.2, yr, 2.3, h, 'RNA\n(100 komponent\nPCA)', RNA, fs=7.2)
    okvir(ax, 0.2, yt, 2.3, h, 'TCR\n(vložitev\nCDR3β)', TCR, fs=7.2)
    trapez(ax, 3.0, yr - 0.1, 2.0, h + 0.2, 'kodirnik\nRNA', RNA)
    trapez(ax, 3.0, yt - 0.1, 2.0, h + 0.2, 'kodirnik\nTCR', TCR)
    puscica(ax, 2.5, yr + h / 2, 3.0, yr + h / 2)
    puscica(ax, 2.5, yt + h / 2, 3.0, yt + h / 2)
    # latentni prostor
    xl, yl, wl, hl = 5.7, 6.25, 3.0, 1.75
    okvir(ax, xl, yl, wl, hl, '', ROB)
    ax.text(xl + wl / 2, yl + hl - 0.3, 'skupni latentni prostor', ha='center', va='center', fontsize=7.5)
    ax.text(xl + wl / 2, yl + hl / 2 - 0.05, r'povprečje $\mu$, $\sigma$ kodirnikov', ha='center', va='center',
            fontsize=7.2)
    ax.text(xl + wl / 2, yl + 0.32, r'vzorec $z$', ha='center', va='center', fontsize=7.5)
    puscica(ax, 5.0, yr + h / 2, xl, yl + hl - 0.45)
    puscica(ax, 5.0, yt + h / 2, xl, yl + 0.45)
    trapez(ax, 9.4, yr - 0.1, 2.0, h + 0.2, 'dekodirnik\nRNA', RNA, ozek_desno=False)
    trapez(ax, 9.4, yt - 0.1, 2.0, h + 0.2, 'dekodirnik\nTCR', TCR, ozek_desno=False)
    puscica(ax, xl + wl, yl + hl - 0.45, 9.4, yr + h / 2)
    puscica(ax, xl + wl, yl + 0.45, 9.4, yt + h / 2)
    okvir(ax, 11.9, yr, 3.9, h, 'rekonstruirana RNA', RNA)
    okvir(ax, 11.9, yt, 3.9, h, 'vektorji TCR: isti klon blizu,\nrazlični kloni narazen', TCR)
    puscica(ax, 11.4, yr + h / 2, 11.9, yr + h / 2)
    puscica(ax, 11.4, yt + h / 2, 11.9, yt + h / 2)
    # pogoji
    xp, yp, wp, hp = 3.3, 9.3, 9.4, 0.75
    okvir(ax, xp, yp, wp, hp, 'pogoji: bolnik, tkivo (kri ali tumor), čas (pred zdravljenjem ali po njem)', SIVA,
          fs=7.2, ls='--')
    for x in (4.0, 10.4):
        puscica(ax, x if x < 5 else x, yp, x, yr + h + 0.1, barva=SIVA, ls='--', lw=0.8)

    # B: napoved
    ax.plot([0.2, 15.8], [4.75, 4.75], color='#e4e3df', lw=0.8)
    crka(ax, 0.1, 4.35, 'B', 'napoved ekspanzije za izpuščenega bolnika')
    yb, hb = 1.15, 2.45
    koraki = [
        (0.2, 2.6, 'celice bolnika\niz krvi pred\nzdravljenjem'),
        (3.35, 2.75, 'kodiranje s\npogojem »kri,\npred zdravljenjem«\n' + r'$\rightarrow$ vzorci $z$'),
        (6.65, 2.75, 'dekodiranje\nvsakega vzorca\ns pogojem »pred«\nin s pogojem »po«'),
        (9.95, 2.75, ''),
        (13.25, 2.55, 'ocena ekspanzije $y_i$:\nvelikost psevdoklona\npo zdravljenju,\ndeljena z velikostjo\npred zdravljenjem'),
    ]
    for i, (x, w, t) in enumerate(koraki):
        okvir(ax, x, yb, w, hb, t, ROB, fs=7.3)
        if i < len(koraki) - 1:
            xn = koraki[i + 1][0]
            puscica(ax, x + w, yb + hb / 2, xn, yb + hb / 2)
    # psevdokloni: ikona
    x4, w4 = koraki[3][0], koraki[3][1]
    ax.text(x4 + w4 / 2, yb + hb - 0.32, 'psevdokloni\ngeneriranih TCR', ha='center', va='top', fontsize=7.3,
            linespacing=1.3)
    rng = np.random.default_rng(3)
    for (cx, cy, n) in [(x4 + 0.7, yb + 0.9, 6), (x4 + 1.45, yb + 0.8, 3), (x4 + 2.1, yb + 0.95, 1)]:
        tocke = rng.normal(0, 0.09, size=(n, 2)) if n > 1 else np.zeros((1, 2))
        ax.scatter(cx + tocke[:, 0], cy + tocke[:, 1], s=6, color=TCR, zorder=3, linewidths=0)
        from matplotlib.patches import Circle
        ax.add_patch(Circle((cx, cy), 0.3, fc='none', ec=SIVA, lw=0.7, ls='--'))
    ax.text(x4 + w4 / 2, yb + 0.22, r'(prag $\tau$)', ha='center', va='center', fontsize=7)
    ax.text(8.0, 0.55, 'celica, katere psevdoklon ima pred zdravljenjem povprečno manj kot dve celici, '
            'dobi oceno 0', ha='center', va='center', fontsize=7, color=ROB)

    os.makedirs(SLIKE, exist_ok=True)
    pot = os.path.join(SLIKE, 'shema_trim.png')
    fig.savefig(pot, facecolor='white', bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(pot)


def shema_scfea(plt):
    from matplotlib.patches import Circle
    s, v = 16.0, 7.4
    fig, ax = platno(plt, s, v)

    # A: modulni graf
    crka(ax, 0.1, 7.05, 'A', 'metabolni moduli in presnovki')
    r = 0.55

    def presnovek(x, y, ime):
        ax.add_patch(Circle((x, y), r, fc='white', ec=ROB, lw=0.9, zorder=3))
        ax.text(x, y, ime, ha='center', va='center', fontsize=6.4, zorder=4)

    def modul(x, y, ime, geni):
        okvir(ax, x - 0.75, y - 0.3, 1.5, 0.6, ime, FLUKS, fs=7.2)
        ax.text(x, y - 0.42, geni, ha='center', va='top', fontsize=6.2, color=ROB, style='italic',
                linespacing=1.2)

    P = {'gluk': (0.8, 4.6), 'pir': (4.05, 4.6), 'lak': (7.3, 5.7), 'acoa': (7.3, 3.3)}
    presnovek(*P['gluk'], 'glukoza')
    presnovek(*P['pir'], 'piruvat')
    presnovek(*P['lak'], 'laktat')
    presnovek(*P['acoa'], 'acetil-\nCoA')
    M = {1: (2.42, 4.6), 2: (5.68, 5.7), 3: (5.68, 3.3)}
    modul(*M[1], 'modul 1', 'HK2, PFKP,\nGAPDH, PKM ...')
    modul(*M[2], 'modul 2', 'LDHA, LDHB')
    modul(*M[3], 'modul 3', 'PDHA1, PDHB,\nDLAT, DLD')
    puscica(ax, P['gluk'][0] + r, 4.6, M[1][0] - 0.75, 4.6)
    puscica(ax, M[1][0] + 0.75, 4.6, P['pir'][0] - r, 4.6)
    puscica(ax, P['pir'][0] + 0.38, 5.0, M[2][0] - 0.75, 5.7)
    puscica(ax, M[2][0] + 0.75, 5.7, P['lak'][0] - r, 5.7)
    puscica(ax, P['pir'][0] + 0.38, 4.2, M[3][0] - 0.75, 3.3)
    puscica(ax, M[3][0] + 0.75, 3.3, P['acoa'][0] - r, 3.3)
    puscica(ax, P['acoa'][0] + r, 3.3, P['acoa'][0] + 1.05, 3.3)
    ax.text(P['acoa'][0] + 1.1, 3.3, '...', ha='left', va='center', fontsize=8)
    ax.text(0.2, 1.65, 'ravnovesje presnovka piruvat:', fontsize=7.3, ha='left', va='center')
    ax.text(0.2, 1.15, r'$f_1 \approx f_2 + f_3$ v vsaki celici;  vsi $f \geq 0$', fontsize=7.6, ha='left',
            va='center')
    ax.text(0.2, 0.6, '168 modulov, povezanih prek 70 presnovkov', fontsize=7.0, ha='left', va='center',
            color=ROB)

    # B: ocena fluksa enega modula
    ax.plot([8.75, 8.75], [0.3, 7.1], color='#e4e3df', lw=0.8)
    crka(ax, 9.0, 7.05, 'B', 'ocena fluksa skozi modul 3')
    okvir(ax, 9.1, 3.6, 2.1, 2.2, 'izražanje\ngenov modula\nv celici $i$\n(PDHA1, PDHB,\nDLAT, DLD)', RNA, fs=7.0)
    # majhna nevronska mreza
    plasti = [(11.95, [4.2, 4.7, 5.2]), (12.75, [3.95, 4.45, 4.95, 5.45]), (13.55, [4.45, 4.95])]
    for (x0, ys0), (x1, ys1) in zip(plasti[:-1], plasti[1:]):
        for y0 in ys0:
            for y1 in ys1:
                ax.plot([x0, x1], [y0, y1], color=SIVA, lw=0.4, zorder=1)
    for x0, ys0 in plasti:
        for y0 in ys0:
            ax.add_patch(Circle((x0, y0), 0.11, fc='white', ec=ROB, lw=0.7, zorder=2))
    ax.text(12.75, 3.35, 'nevronska mreža\nmodula 3', ha='center', va='top', fontsize=7.0, linespacing=1.25)
    puscica(ax, 11.2, 4.7, 11.8, 4.7)
    puscica(ax, 13.7, 4.7, 14.3, 4.7)
    okvir(ax, 14.3, 4.2, 1.55, 1.0, r'fluks $f_{i3}$', FLUKS, fs=7.4)
    ax.text(12.45, 1.65, 'mreže vseh modulov se učijo skupaj, brez meritev fluksov:', ha='center', va='center',
            fontsize=7.0)
    ax.text(12.45, 1.25, 'ravnovesje presnovkov, nenegativni fluksi in povezanost', ha='center', va='center',
            fontsize=7.0)
    ax.text(12.45, 0.85, 'fluksov z izražanjem genov modulov', ha='center', va='center', fontsize=7.0)

    pot = os.path.join(SLIKE, 'shema_scfea.png')
    fig.savefig(pot, facecolor='white', bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(pot)


def main():
    plt = pripravi_plt()
    shema_trim(plt)
    shema_scfea(plt)


if __name__ == '__main__':
    main()
