"""Slika 10: ekspanzija glede na velikost klona pred zdravljenjem in globino vzorčenja.

Bere izhode zvezka colab_preverba_klonov.ipynb v diploma/rezultati/preverba_klonov/:
    velikosti_klonov.csv      delež ekspandiranih klonov po velikosti klona v krvi pred zdravljenjem (bolnika 23, 24)
    globina_po_bolnikih.csv   delež ekspandiranih celic in razmerje števila celic krvi po/pred zdravljenju (vsi bolniki)

Zagon iz korena repozitorija:
    python diploma/skripte/slika_velikost_klona.py

Izhod:
    diploma/slike/velikost_klona.png
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

KOREN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
VHOD = os.path.join(KOREN, 'diploma', 'rezultati', 'preverba_klonov')
SLIKE = os.path.join(KOREN, 'diploma', 'slike')
VELIKOSTI = ['1', '2', '3-5', '6-10', '>10']
OZNAKE_VEL = {'1': '1', '2': '2', '3-5': '3–5', '6-10': '6–10', '>10': '>10'}
STOLPEC = '#52514e'
POUDAREK = '#2a78d6'
SIVA = '#a3a29b'
MREZA = '#e4e3df'


def vejica(x, _=None, dec=1):
    return f'{x:.{dec}f}'.replace('.', ',')


def pripravi_plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
        'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e',
        'ytick.color': '#52514e', 'axes.linewidth': 0.8, 'savefig.dpi': 300,
        'axes.formatter.use_locale': False,
    })
    return plt


def main():
    from matplotlib.ticker import FuncFormatter
    plt = pripravi_plt()
    vel = pd.read_csv(os.path.join(VHOD, 'velikosti_klonov.csv'), dtype={'velikost_krv_pred': str})
    glob_ = pd.read_csv(os.path.join(VHOD, 'globina_po_bolnikih.csv'))
    glob_ = glob_[(glob_.krv_pred > 0) & (glob_.krv_po > 0)].copy()
    rho, p = spearmanr(glob_.razmerje_po_pred, glob_.delez_eksp_vsi_bolniki)
    rho_d, p_d = spearmanr(glob_.razmerje_po_pred, glob_.delez_eksp_delezi)

    fig, osi = plt.subplots(1, 3, figsize=(16 / 2.54, 6.4 / 2.54),
                            gridspec_kw={'width_ratios': [1, 1, 1.15], 'wspace': 0.42})
    for ax, crka, bolnik in zip(osi[:2], 'AB', (23, 24)):
        sub = vel[vel.bolnik == bolnik].set_index('velikost_krv_pred')
        x = np.arange(len(VELIKOSTI))
        for i, v in enumerate(VELIKOSTI):
            if v not in sub.index:
                ax.text(i, 0.03, '–', ha='center', va='bottom', fontsize=7.5, color=SIVA)
                continue
            r = sub.loc[v]
            d = float(r.delez_eksp_klonov)
            ax.bar(i, d, width=0.62, color=STOLPEC, zorder=2)
            ax.text(i, d + 0.03, f'{int(r.n_eksp_klonov)}/{int(r.n_klonov)}', ha='center', va='bottom',
                    fontsize=7, color='#0b0b0b')
        ax.set_xticks(x)
        ax.set_xticklabels([OZNAKE_VEL[v] for v in VELIKOSTI], fontsize=7)
        ax.set_ylim(0, 1.12)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.yaxis.set_major_formatter(FuncFormatter(vejica))
        ax.grid(axis='y', color=MREZA, lw=0.6, zorder=0)
        ax.tick_params(axis='x', length=0)
        ax.set_xlabel('celic klona v krvi pred zdr.')
        ax.text(-0.02, 1.06, crka, transform=ax.transAxes, fontweight='bold', fontsize=10, va='bottom', ha='right')
        ax.text(0.0, 1.06, f'bolnik {bolnik}', transform=ax.transAxes, fontsize=8.5, va='bottom')
    osi[0].set_ylabel('delež ekspandiranih klonov')

    ax = osi[2]
    druge = glob_[~glob_.bolnik.isin([23, 24])]
    ax.scatter(druge.razmerje_po_pred, druge.delez_eksp_vsi_bolniki, s=16, color=SIVA, edgecolor='white',
               linewidth=0.5, zorder=2)
    izp = glob_[glob_.bolnik.isin([23, 24])]
    ax.scatter(izp.razmerje_po_pred, izp.delez_eksp_vsi_bolniki, s=22, color=POUDAREK, edgecolor='white',
               linewidth=0.5, zorder=3)
    for _, r in izp.iterrows():
        ax.annotate(f'{int(r.bolnik)}', (r.razmerje_po_pred, r.delez_eksp_vsi_bolniki), xytext=(4, -3),
                    textcoords='offset points', fontsize=7.5, ha='left', va='top', color='#0b0b0b')
    ax.set_xlim(0.3, 1.6)
    ax.set_ylim(-0.01, 0.34)
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_xticks([0.4, 0.8, 1.2, 1.6])
    ax.xaxis.set_major_formatter(FuncFormatter(vejica))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _=None: vejica(v, dec=1)))
    ax.grid(color=MREZA, lw=0.6, zorder=0)
    ax.set_xlabel('celice v krvi po/pred zdr.')
    ax.set_ylabel('delež ekspandiranih celic')
    ax.text(0.03, 0.97, f'ρ = {vejica(rho, dec=2)}', transform=ax.transAxes, fontsize=7.5, va='top')
    ax.text(-0.02, 1.06, 'C', transform=ax.transAxes, fontweight='bold', fontsize=10, va='bottom', ha='right')
    ax.text(0.0, 1.06, f'{len(glob_)} bolnikov', transform=ax.transAxes, fontsize=8.5, va='bottom')

    os.makedirs(SLIKE, exist_ok=True)
    pot = os.path.join(SLIKE, 'velikost_klona.png')
    fig.savefig(pot, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(pot)
    print(json.dumps({'n_bolnikov': int(len(glob_)), 'spearman_stevila': [round(rho, 3), round(p, 4)],
                      'spearman_delezi': [round(rho_d, 3), round(p_d, 4)]}))


if __name__ == '__main__':
    main()
