"""Rezultati napovedi klonske ekspanzije iz zapisov zagonov.

Bere results/logs/logs/*.json (resnične oznake x_true in ocene y_score za vsako celico)
in izračuna vse številke za poglavje Rezultati: AUC z DeLongovim 95-% intervalom zaupanja,
povzetke po semenih, parne DeLongove teste med različicami in med semeni, metrike pri
pragovih ter primerjavo z izhodiščnimi klasifikatorji. Nariše sliki z krivuljami ROC in z AUC po semenih.

Zagon iz korena repozitorija:
    python diploma/skripte/rezultati_napovedi.py

Izhodi:
    diploma/rezultati/auc_zagoni.csv        en zagon na vrstico
    diploma/rezultati/auc_povzetek.csv      različica x bolnik (povprečje, razpon po semenih)
    diploma/rezultati/delong_razlicice.csv  parni DeLongov test med različicami (isti bolnik, isto seme)
    diploma/rezultati/delong_semena.csv     parni DeLongov test med semeni (ista različica, isti bolnik)
    diploma/rezultati/pragovi.csv           metrike pri treh pragovih (iz zapisov)
    diploma/rezultati/povzetek.json         številke, ki jih navaja besedilo
    diploma/rezultati/tabele.md             tabele 4–6 za poglavje Rezultati
    diploma/slike/roc_izpusceni_bolniki.png  krivulje ROC (slika 6)
    diploma/slike/auc_semena_izhodisca.png   AUC po semenih z intervali zaupanja in izhodišči (slika 7)

Ničesar zunaj mape diploma/ ne spreminja.
"""
from __future__ import annotations

import glob
import json
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr
from sklearn.metrics import roc_auc_score, roc_curve

KOREN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
LOGI = os.path.join(KOREN, 'results', 'logs', 'logs')
IZHOD = os.path.join(KOREN, 'diploma', 'rezultati')
SLIKE = os.path.join(KOREN, 'diploma', 'slike')

RAZLICICE = ['rna_tcr', 'flux_tcr', 'rna_flux_tcr']
IMENA = {'rna_tcr': 'RNA + TCR', 'flux_tcr': 'fluksi + TCR', 'rna_flux_tcr': 'RNA + fluksi + TCR'}
# kategorične barve (prve tri reže preverjene palete; ločljive tudi pri barvni slepoti) + slog črte za črno-beli tisk
BARVE = {'rna_tcr': '#2a78d6', 'flux_tcr': '#eb6834', 'rna_flux_tcr': '#1baf7a'}
CRTE = {'rna_tcr': '-', 'flux_tcr': '--', 'rna_flux_tcr': ':'}
OZNAKE = {'rna_tcr': 'o', 'flux_tcr': 's', 'rna_flux_tcr': '^'}
IZHODISCA = ['knn_rna', 'NN_rna', 'knn_flux', 'NN_flux', 'knn_tcr', 'NN_tcr', 'knn_both', 'NN_both']
MEJA_TAU = 10 * 0.95 ** 50          # umerjanje praga se je ustavilo pri 50 ponovitvah
Z95 = norm.ppf(0.975)


# --------------------------------------------------------------------------- DeLong
def delong(y: np.ndarray, ocene: np.ndarray):
    """Hitri DeLongov algoritem (Sun in Xu, 2014) za k korelirano ocen na istih celicah.

    y: 0/1 dolžine n; ocene: k x n. Vrne (AUC dolžine k, kovariančna matrika k x k).
    Vezi dobijo povprečni rang, kar ustreza AUC iz sklearn.
    """
    ocene = np.atleast_2d(np.asarray(ocene, dtype=float))
    y = np.asarray(y).astype(bool)
    poz, neg = ocene[:, y], ocene[:, ~y]
    m, n = poz.shape[1], neg.shape[1]
    tx = np.vstack([rankdata(r) for r in poz])
    ty = np.vstack([rankdata(r) for r in neg])
    tz = np.vstack([rankdata(r) for r in np.hstack([poz, neg])])
    auc = tz[:, :m].sum(axis=1) / (m * n) - (m + 1) / (2 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1 - (tz[:, m:] - ty) / m
    s = np.atleast_2d(np.cov(v01)) / m + np.atleast_2d(np.cov(v10)) / n
    return auc, s


def delong_iz(y, ocena):
    auc, s = delong(y, ocena)
    se = float(np.sqrt(s[0, 0]))
    return float(auc[0]), se, max(0.0, auc[0] - Z95 * se), min(1.0, auc[0] + Z95 * se)


def delong_par(y, a, b):
    auc, s = delong(y, np.vstack([a, b]))
    var = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    z = (auc[0] - auc[1]) / np.sqrt(var) if var > 0 else 0.0
    return float(auc[0]), float(auc[1]), float(z), float(2 * norm.sf(abs(z)))


# --------------------------------------------------------------------------- branje
def preberi_zapise():
    zapisi = []
    for pot in sorted(glob.glob(os.path.join(LOGI, '*.json'))):
        with open(pot, encoding='utf-8') as f:
            d = json.load(f)
        d['_datoteka'] = os.path.basename(pot)
        zapisi.append(d)
    # Zagon flux_tcr, bolnik 23, seme 0 je zapisan dvakrat (ponovitev za preverjanje ponovljivosti).
    # V glavni analizi uporabimo PRVI zagon (najzgodnejši časovni žig), drugega samo za ponovljivost.
    zapisi.sort(key=lambda d: (d['variant'], d['heldout_patient'], d['seed'], d['timestamp']))
    videno = set()
    for d in zapisi:
        kljuc = (d['variant'], d['heldout_patient'], d['seed'])
        d['_glavni'] = kljuc not in videno
        videno.add(kljuc)
    return zapisi


# --------------------------------------------------------------------------- izračuni
def tabela_zagonov(zapisi):
    vrstice = []
    for d in zapisi:
        y = np.asarray(d['x_true'], dtype=int)
        s = np.asarray(d['y_score'], dtype=float)
        auc, se, lo, hi = delong_iz(y, s)
        assert abs(roc_auc_score(y, s) - auc) < 1e-9, d['_datoteka']
        assert abs(auc - d['roc_auc']) < 6e-5, (d['_datoteka'], auc, d['roc_auc'])
        izh = {k: v for k, v in (d.get('baseline_results') or {}).items()}
        najb_ime = max(izh, key=izh.get) if izh else None
        vrstica = {
            'razlicica': d['variant'], 'bolnik': d['heldout_patient'], 'seme': d['seed'],
            'casovni_zig': d['timestamp'], 'glavni': d['_glavni'], 'datoteka': d['_datoteka'],
            'n_celic': len(y), 'n_ekspandiranih': int(y.sum()), 'delez_ekspandiranih': y.mean(),
            'auc_zapis': d['roc_auc'], 'auc': auc, 'se': se, 'ip_spodaj': lo, 'ip_zgoraj': hi,
            'tau': d['thresh_fitted'], 'tau_na_meji': d['thresh_fitted'] < MEJA_TAU + 1e-6,
            'delez_y0': float((s == 0).mean()),
            'delez_eksp_med_y0': float(y[s == 0].mean()) if (s == 0).any() else np.nan,
            'delez_eksp_med_y_poz': float(y[s > 0].mean()) if (s > 0).any() else np.nan,
            'najboljse_izhodisce': najb_ime, 'auc_najboljse_izhodisce': izh.get(najb_ime, np.nan),
        }
        for k in IZHODISCA:
            vrstica[k] = izh.get(k, np.nan)
        vrstica['ip_nad_izhodisci'] = lo > vrstica['auc_najboljse_izhodisce']
        vrstice.append(vrstica)
    return pd.DataFrame(vrstice)


def tabela_povzetka(df):
    g = df[df.glavni].groupby(['razlicica', 'bolnik'])
    p = g['auc'].agg(n_semen='count', povprecje='mean', sd='std', najmanj='min', najvec='max').reset_index()
    p['izhodisce_najvec'] = g['auc_najboljse_izhodisce'].max().values
    p['zagonov_ip_nad_izhodisci'] = g['ip_nad_izhodisci'].sum().values
    p['zagonov_tau_na_meji'] = g['tau_na_meji'].sum().values
    p['razlicica'] = pd.Categorical(p['razlicica'], RAZLICICE, ordered=True)
    return p.sort_values(['bolnik', 'razlicica']).reset_index(drop=True)


def testi_razlicic(zapisi):
    glavni = {(d['variant'], d['heldout_patient'], d['seed']): d for d in zapisi if d['_glavni']}
    pari = [('rna_flux_tcr', 'rna_tcr'), ('flux_tcr', 'rna_tcr'), ('rna_flux_tcr', 'flux_tcr')]
    vrstice = []
    for bolnik in sorted({k[1] for k in glavni}):
        for seme in sorted({k[2] for k in glavni}):
            for a, b in pari:
                if (a, bolnik, seme) not in glavni or (b, bolnik, seme) not in glavni:
                    continue
                da, db = glavni[(a, bolnik, seme)], glavni[(b, bolnik, seme)]
                assert da['x_true'] == db['x_true'], 'Oznake se ne ujemajo – celice niso iste.'
                y = np.asarray(da['x_true'])
                auc_a, auc_b, z, p = delong_par(y, da['y_score'], db['y_score'])
                vrstice.append({'bolnik': bolnik, 'seme': seme, 'a': a, 'b': b, 'auc_a': auc_a,
                                'auc_b': auc_b, 'razlika': auc_a - auc_b, 'z': z, 'p': p})
    return pd.DataFrame(vrstice)


def testi_semen(zapisi):
    vrstice = []
    skupine = {}
    for d in zapisi:
        skupine.setdefault((d['variant'], d['heldout_patient']), []).append(d)
    for (razl, bolnik), ds in sorted(skupine.items()):
        for da, db in combinations(ds, 2):
            y = np.asarray(da['x_true'])
            auc_a, auc_b, z, p = delong_par(y, da['y_score'], db['y_score'])
            vrstice.append({'razlicica': razl, 'bolnik': bolnik,
                            'oba_glavna': bool(da['_glavni'] and db['_glavni']),
                            'zagon_a': f"seme {da['seed']} ({da['timestamp']})",
                            'zagon_b': f"seme {db['seed']} ({db['timestamp']})",
                            'isto_seme': da['seed'] == db['seed'],
                            'auc_a': auc_a, 'auc_b': auc_b, 'razlika': auc_a - auc_b, 'z': z, 'p': p,
                            'spearman_ocen': float(spearmanr(da['y_score'], db['y_score'])[0])})
    return pd.DataFrame(vrstice)


def tabela_pragov(zapisi):
    vrstice = []
    for d in zapisi:
        if not d['_glavni']:
            continue
        for ime, t in (d.get('thresholds') or {}).items():
            c = t['confusion']
            vrstice.append({'razlicica': d['variant'], 'bolnik': d['heldout_patient'], 'seme': d['seed'],
                            'prag': ime, 'vrednost_praga': t['thresh'],
                            'n_napovedanih_eksp': t['n_pred_expanded'], **c,
                            'priklic': t['EXPANDED']['recall'], 'preciznost': t['EXPANDED']['precision'],
                            'f1': t['EXPANDED']['f1'], 'specificnost': t['specificity']})
    return pd.DataFrame(vrstice)


# --------------------------------------------------------------------------- tabele za besedilo
def st(x, dec=3):
    return f'{x:.{dec}f}'.replace('.', ',')


def p_besedilo(p):
    return '<\u00a00,001' if p < 0.001 else st(p, 3 if p < 0.1 else 2)


def tabele_md(df, povz, raz):
    glavni = df[df.glavni]
    vrst = ['Table: Tabela 4: AUC napovedi klonske ekspanzije modela TRIM pri izpuščenih bolnikih 23 in 24. '
            'V oklepaju je 95-% interval zaupanja po DeLongu; »–« pomeni, da zagona nismo izvedli. '
            'Povprečje je izračunano po razpoložljivih semenih.', '',
            '| Različica | Seme | Bolnik 23 | Bolnik 24 |',
            '|----------------------|-------------|--------------------------|--------------------------|']
    for razl in RAZLICICE:
        for i, seme in enumerate((0, 1, 2)):
            celice = []
            for b in (23, 24):
                r = glavni[(glavni.razlicica == razl) & (glavni.bolnik == b) & (glavni.seme == seme)]
                celice.append('–' if r.empty else f'{st(r.auc.iloc[0])} ({st(r.ip_spodaj.iloc[0])}–{st(r.ip_zgoraj.iloc[0])})')
            vrst.append(f"| {IMENA[razl] if i == 0 else ''} | {seme} | " + ' | '.join(celice) + ' |')
        pv = [st(povz[(povz.razlicica == razl) & (povz.bolnik == b)].povprecje.iloc[0]) for b in (23, 24)]
        vrst.append(f'|  | povprečje | {pv[0]} | {pv[1]} |')

    vrst += ['', 'Table: Tabela 5: AUC izhodiščnih klasifikatorjev. Pri vhodih, ki smo jih uporabili v več zagonih, '
             'je naveden razpon vrednosti.', '', '| Klasifikator | Vhod | Bolnik 23 | Bolnik 24 |', '|---|---|---|---|']
    vhodi = [('rna', 'RNA', ['rna_tcr', 'rna_flux_tcr']), ('tcr', 'TCR', RAZLICICE),
             ('both', 'RNA in TCR', ['rna_tcr', 'rna_flux_tcr']), ('flux', 'fluksi', ['flux_tcr']),
             ('both', 'fluksi in TCR', ['flux_tcr'])]
    for klas, ime_k in (('knn', 'kNN'), ('NN', 'NN')):
        for kljuc, ime_v, razl in vhodi:
            celice = []
            for b in (23, 24):
                v = glavni[(glavni.bolnik == b) & glavni.razlicica.isin(razl)][f'{klas}_{kljuc}'].dropna()
                lo, hi = round(v.min(), 3), round(v.max(), 3)
                celice.append(st(lo) if lo == hi else f'{st(lo)}–{st(hi)}')
            vrst.append(f'| {ime_k} | {ime_v} | ' + ' | '.join(celice) + ' |')

    vrst += ['', 'Table: Tabela 6: Razlika v AUC (ΔAUC) različic RNA + fluksi + TCR in fluksi + TCR do različice '
             'RNA + TCR na istih celicah pri istem bolniku in semenu, s parnim DeLongovim testom. Pozitivna razlika '
             'pomeni višjo AUC različice s fluksi; pri bolniku 24 s semenom 2 zagona RNA + TCR nismo izvedli.', '',
             '| Bolnik | Seme | RNA + fluksi + TCR: ΔAUC | *p* | fluksi + TCR: ΔAUC | *p* |',
             '|----------|--------|------------------------|------------|------------------------|------------|']
    for (b, seme), g in raz.groupby(['bolnik', 'seme']):
        a = g[(g.a == 'rna_flux_tcr') & (g.b == 'rna_tcr')]
        f = g[(g.a == 'flux_tcr') & (g.b == 'rna_tcr')]
        if a.empty and f.empty:
            continue
        def cel(x):
            return ('–', '–') if x.empty else (st(x.razlika.iloc[0]).replace('-', '−'), p_besedilo(x.p.iloc[0]))
        ca, cf = cel(a), cel(f)
        vrst.append(f'| {b} | {seme} | {ca[0]} | {ca[1]} | {cf[0]} | {cf[1]} |')
    return '\n'.join(vrst) + '\n'


# --------------------------------------------------------------------------- slike
def nastavi_matplotlib():
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


def vejica(x, _=None, dec=1):
    return f'{x:.{dec}f}'.replace('.', ',')


def slika_roc(zapisi, plt, pot):
    from matplotlib.ticker import FuncFormatter
    from matplotlib.lines import Line2D
    bolniki = sorted({d['heldout_patient'] for d in zapisi})
    fig, osi = plt.subplots(1, len(bolniki), figsize=(16 / 2.54, 7.6 / 2.54), sharey=True)
    for ax, crka, bolnik in zip(osi, 'AB', bolniki):
        ax.plot([0, 1], [0, 1], color='#a3a29b', lw=0.8, ls=(0, (2, 2)), zorder=1)
        n = n1 = 0
        for razl in RAZLICICE:
            for d in zapisi:
                if d['variant'] != razl or d['heldout_patient'] != bolnik or not d['_glavni']:
                    continue
                fpr, tpr, _ = roc_curve(d['x_true'], d['y_score'])
                ax.plot(fpr, tpr, color=BARVE[razl], ls=CRTE[razl], lw=1.3, alpha=0.9, zorder=2)
                n, n1 = len(d['x_true']), int(np.sum(d['x_true']))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.01); ax.set_aspect('equal')
        ax.set_xticks(np.linspace(0, 1, 6)); ax.set_yticks(np.linspace(0, 1, 6))
        ax.xaxis.set_major_formatter(FuncFormatter(vejica)); ax.yaxis.set_major_formatter(FuncFormatter(vejica))
        ax.set_xlabel('delež lažno pozitivnih')
        ax.text(-0.02, 1.08, crka, transform=ax.transAxes, fontweight='bold', fontsize=10, va='bottom', ha='right')
        ax.text(0.0, 1.08, f'bolnik {bolnik} ({n} celic, {n1} ekspandiranih)', transform=ax.transAxes,
                fontsize=8.5, va='bottom', color='#0b0b0b')
    osi[0].set_ylabel('delež pravilno pozitivnih')
    rocice = [Line2D([], [], color=BARVE[r], ls=CRTE[r], lw=1.5, label=IMENA[r]) for r in RAZLICICE]
    rocice.append(Line2D([], [], color='#a3a29b', ls=(0, (2, 2)), lw=0.8, label='naključno ugibanje'))
    fig.legend(handles=rocice, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01),
               handlelength=2.6, columnspacing=1.4)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(pot, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def slika_auc(df, plt, pot):
    from matplotlib.ticker import FuncFormatter
    from matplotlib.lines import Line2D
    glavni = df[df.glavni]
    bolniki = sorted(glavni.bolnik.unique())
    fig, osi = plt.subplots(1, len(bolniki), figsize=(16 / 2.54, 6.8 / 2.54), sharey=True)
    for ax, crka, bolnik in zip(osi, 'AB', bolniki):
        ax.axhline(0.5, color='#a3a29b', lw=0.8, ls=(0, (2, 2)), zorder=1)
        for i, razl in enumerate(RAZLICICE):
            sub = glavni[(glavni.bolnik == bolnik) & (glavni.razlicica == razl)].sort_values('seme')
            # najboljši izhodiščni klasifikator (največja AUC med šestimi) za to različico: siva vodoravna črta
            for _, r in sub.iterrows():
                x = i + (r.seme - 1) * 0.2
                ax.hlines(r.auc_najboljse_izhodisce, x - 0.08, x + 0.08, color='#52514e', lw=1.6, zorder=2)
                ax.errorbar(x, r.auc, yerr=[[r.auc - r.ip_spodaj], [r.ip_zgoraj - r.auc]], fmt=OZNAKE[razl],
                            color=BARVE[razl], ms=5.5, mec='white', mew=0.8, elinewidth=1.2, capsize=0, zorder=3)
                ax.text(x, 0.385, f'{int(r.seme)}', ha='center', va='center', fontsize=7, color='#52514e')
        ax.set_xticks(range(len(RAZLICICE)))
        ax.set_xticklabels([IMENA[r].replace(' + fluksi + ', ' +\nfluksi + ') for r in RAZLICICE])
        ax.set_xlim(-0.5, len(RAZLICICE) - 0.5)
        ax.set_ylim(0.36, 1.0)
        ax.set_yticks(np.arange(0.4, 1.01, 0.1))
        ax.yaxis.set_major_formatter(FuncFormatter(vejica))
        ax.grid(axis='y', color='#e4e3df', lw=0.6, zorder=0)
        ax.tick_params(axis='x', length=0)
        ax.text(-0.02, 1.04, crka, transform=ax.transAxes, fontweight='bold', fontsize=10, va='bottom', ha='right')
        ax.text(0.0, 1.04, f'bolnik {bolnik}', transform=ax.transAxes, fontsize=8.5, va='bottom')
    osi[0].set_ylabel('AUC')
    rocice = [Line2D([], [], marker=OZNAKE[r], color=BARVE[r], ls='', ms=5.5, mec='white', label=IMENA[r])
              for r in RAZLICICE]
    rocice.append(Line2D([], [], color='#52514e', lw=1.6, label='najboljši izhodiščni klasifikator'))
    rocice.append(Line2D([], [], color='#a3a29b', lw=0.8, ls=(0, (2, 2)), label='naključno ugibanje'))
    rocice = [rocice[i] for i in (0, 3, 1, 4, 2)]   # legenda se bere po vrsticah
    fig.legend(handles=rocice, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.03),
               columnspacing=1.2, handlelength=1.8)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(pot, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# --------------------------------------------------------------------------- glavno
def main():
    os.makedirs(IZHOD, exist_ok=True)
    os.makedirs(SLIKE, exist_ok=True)
    zapisi = preberi_zapise()
    df = tabela_zagonov(zapisi)
    povz = tabela_povzetka(df)
    raz = testi_razlicic(zapisi)
    sem = testi_semen(zapisi)
    prag = tabela_pragov(zapisi)

    df.to_csv(os.path.join(IZHOD, 'auc_zagoni.csv'), index=False, float_format='%.4f')
    povz.to_csv(os.path.join(IZHOD, 'auc_povzetek.csv'), index=False, float_format='%.4f')
    raz.to_csv(os.path.join(IZHOD, 'delong_razlicice.csv'), index=False, float_format='%.4g')
    sem.to_csv(os.path.join(IZHOD, 'delong_semena.csv'), index=False, float_format='%.4g')
    prag.to_csv(os.path.join(IZHOD, 'pragovi.csv'), index=False, float_format='%.4f')

    glavni = df[df.glavni]
    rep = sem[sem.isto_seme].iloc[0].to_dict() if sem.isto_seme.any() else None
    razl_sem = sem[~sem.isto_seme]
    povzetek = {
        'zagonov_vseh': int(len(df)), 'zagonov_glavnih': int(len(glavni)),
        'tau_na_meji': int(glavni.tau_na_meji.sum()), 'meja_tau': MEJA_TAU,
        'tau_na_meji_zagoni': glavni.loc[glavni.tau_na_meji, ['razlicica', 'bolnik', 'seme']].to_dict('records'),
        'ponovljeni_zagon': rep,
        'pari_semen_p_pod_0_05': int((razl_sem.p < 0.05).sum()), 'pari_semen_vseh': int(len(razl_sem)),
        'glavni_pari_semen_p_pod_0_05': int((razl_sem[razl_sem.oba_glavna].p < 0.05).sum()),
        'glavni_pari_semen_vseh': int(razl_sem.oba_glavna.sum()),
        'spearman_med_semeni': {r: [float(razl_sem[razl_sem.oba_glavna & (razl_sem.razlicica == r)].spearman_ocen.min()),
                                    float(razl_sem[razl_sem.oba_glavna & (razl_sem.razlicica == r)].spearman_ocen.max())]
                                for r in RAZLICICE},
        'delez_y0': {f'{r}_{b}': [float(g.delez_y0.min()), float(g.delez_y0.max())]
                     for (r, b), g in glavni.groupby(['razlicica', 'bolnik'])},
        'razlicice_p_pod_0_05': int((raz.p < 0.05).sum()), 'razlicice_vseh': int(len(raz)),
    }
    with open(os.path.join(IZHOD, 'povzetek.json'), 'w', encoding='utf-8') as f:
        json.dump(povzetek, f, ensure_ascii=False, indent=2, default=float)

    with open(os.path.join(IZHOD, 'tabele.md'), 'w', encoding='utf-8') as f:
        f.write(tabele_md(df, povz, raz))

    plt = nastavi_matplotlib()
    slika_roc(zapisi, plt, os.path.join(SLIKE, 'roc_izpusceni_bolniki.png'))
    slika_auc(df, plt, os.path.join(SLIKE, 'auc_semena_izhodisca.png'))

    pd.set_option('display.width', 200)
    print(glavni[['razlicica', 'bolnik', 'seme', 'auc', 'ip_spodaj', 'ip_zgoraj', 'tau', 'delez_y0',
                  'najboljse_izhodisce', 'auc_najboljse_izhodisce', 'ip_nad_izhodisci']].round(3).to_string(index=False))
    print(); print(povz.round(3).to_string(index=False))
    print(); print(raz.round(4).to_string(index=False))
    print(); print(sem.round(4).to_string(index=False))
    print(); print(json.dumps(povzetek, ensure_ascii=False, indent=1, default=float))


if __name__ == '__main__':
    main()
