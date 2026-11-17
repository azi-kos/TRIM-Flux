"""Intervali zaupanja po klonih za vse zagone modela TRIM in primerjava na celicah izhodiščnih klasifikatorjev.

Celice istega klona imajo isto oznako ekspanzije in podobne ocene, zato DeLongovi intervali na celicah
negotovost podcenijo. Skripta za vsak zagon izračuna 95-% interval bootstrap po klonih (2000 ponovitev,
v vsaki s ponavljanjem izberemo klone in vzamemo vse njihove celice) ter parne primerjave
(različici pri istem semenu, semeni pri isti različici) z istim bootstrapom.

Izhodiščna klasifikatorja sta bila testirana samo na celicah izpuščenega bolnika iz krvi pred zdravljenjem,
katerih zaporedja CDR3β ni med celicami iz krvi ostalih bolnikov (notebooks/04_expansion_eval.ipynb,
baseline_expansion_prediction_single). Skripta zato AUC modela TRIM izračuna tudi na tej podmnožici,
da ga lahko primerjamo z izhodiščnima klasifikatorjema na istih celicah.

Potrebuje data_labels.pkl in df_all_tcrs.pkl (Drive: MyDrive/Diploma/data/processed; za branje pip install pyarrow)
in zapise zagonov v results/logs/logs.

Zagon iz korena repozitorija:
    PREVERBA_PROC=<mapa s pickli> python diploma/skripte/intervali_po_klonih.py

Izhodi (diploma/rezultati/intervali_po_klonih/):
    zagoni.csv            AUC vsakega zagona z intervaloma DeLong in bootstrap po klonih, na vseh celicah
                          in na celicah izhodiščnih klasifikatorjev, ter najboljši izhodiščni klasifikator
    razlicice.csv         ΔAUC različic do RNA + TCR pri istem bolniku in semenu (parni bootstrap po klonih)
    semena.csv            ΔAUC med semeni iste različice (parni bootstrap po klonih)
    povzetek.json         številke za besedilo
"""
from __future__ import annotations

import glob
import json
import os
import pickle
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

KOREN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
PROC = os.environ.get('PREVERBA_PROC', '/content/drive/MyDrive/Diploma/data/processed')
LOGI = os.environ.get('PREVERBA_LOGI', os.path.join(KOREN, 'results', 'logs', 'logs'))
IZHOD = os.path.join(KOREN, 'diploma', 'rezultati', 'intervali_po_klonih')
N_BOOT, SEME = 2000, 0
Z95 = norm.ppf(0.975)
RAZLICICE = ['rna_tcr', 'flux_tcr', 'rna_flux_tcr']


def auc_stolpci(x, M):
    x = np.asarray(x, bool); M = np.asarray(M, float).reshape(len(x), -1)
    m = int(x.sum()); n = len(x) - m
    if m == 0 or n == 0:
        return np.full(M.shape[1], np.nan)
    r = rankdata(M, axis=0)
    return (r[x].sum(0) - m * (m + 1) / 2) / (m * n)


def delong(x, M):
    x = np.asarray(x, bool); M = np.asarray(M, float).reshape(len(x), -1).T
    poz, neg = M[:, x], M[:, ~x]; m, n = poz.shape[1], neg.shape[1]
    tx = rankdata(poz, axis=1); ty = rankdata(neg, axis=1); tz = rankdata(np.hstack([poz, neg]), axis=1)
    auc = tz[:, :m].sum(1) / (m * n) - (m + 1) / (2 * n)
    v01 = (tz[:, :m] - tx) / n; v10 = 1 - (tz[:, m:] - ty) / m
    S = np.atleast_2d(np.cov(v01)) / m + np.atleast_2d(np.cov(v10)) / n
    return auc, S


def boot_auc(x, M, kloni):
    """AUC stolpcev M v N_BOOT vzorcih klonov s ponavljanjem (ista izbira klonov za vse stolpce)."""
    x = np.asarray(x, bool); M = np.asarray(M, float).reshape(len(x), -1)
    _, inv = np.unique(kloni, return_inverse=True); inv = inv.ravel()
    red = np.argsort(inv, kind='stable'); dolz = np.bincount(inv)
    zac = np.concatenate([[0], np.cumsum(dolz)[:-1]]); K = len(dolz)
    rng = np.random.default_rng(SEME)
    out = np.empty((N_BOOT, M.shape[1]))
    for b in range(N_BOOT):
        s = rng.integers(0, K, K); L = dolz[s]
        idx = red[np.repeat(zac[s] - np.cumsum(L) + L, L) + np.arange(L.sum())]
        out[b] = auc_stolpci(x[idx], M[idx])
    return out


def interval(v):
    v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [2.5, 97.5])) if len(v) else (np.nan, np.nan)


def p_boot(d):
    d = d[np.isfinite(d)]
    return float(min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))) if len(d) else np.nan


def nalozi_zapise():
    zapisi = []
    for pot in sorted(glob.glob(os.path.join(LOGI, '*.json'))):
        r = json.load(open(pot, encoding='utf-8'))
        if r.get('seed') is None:
            continue
        r['_datoteka'] = os.path.basename(pot)
        zapisi.append(r)
    # pri podvojenem zagonu je glavni prvi po času (kot rezultati_napovedi.py)
    zapisi.sort(key=lambda r: (r['variant'], r['heldout_patient'], r['seed'], str(r['timestamp'])))
    videni = set()
    for r in zapisi:
        k = (r['variant'], r['heldout_patient'], r['seed'])
        r['_glavni'] = k not in videni
        videni.add(k)
    return zapisi


def main():
    with open(os.path.join(PROC, 'data_labels.pkl'), 'rb') as f:
        data_labels = pickle.load(f)
    with open(os.path.join(PROC, 'df_all_tcrs.pkl'), 'rb') as f:
        df_all_tcrs = pickle.load(f)
    c_t = data_labels.columns.get_loc('Tissue'); c_s = data_labels.columns.get_loc('Treatment Stage')
    c_p = data_labels.columns.get_loc('Patient'); c_c = data_labels.columns.get_loc('CDR3(Beta1)')
    xl = data_labels.values
    xl = xl[xl[:, c_c] != -1]
    stevila = df_all_tcrs.fillna(0).values

    zapisi = nalozi_zapise()
    os.makedirs(IZHOD, exist_ok=True)
    v_zag, v_razl, v_sem, povz = [], [], [], {}
    for pid in sorted({r['heldout_patient'] for r in zapisi}):
        ho = xl[xl[:, c_p] == pid]
        bb = (ho[:, c_t] == 0) & (ho[:, c_s] == 0)
        tcr = ho[bb, c_c].astype(int)
        x = stevila[tcr, 0] < stevila[tcr, 1]
        pre = stevila[tcr, 0]
        # celice izhodiščnih klasifikatorjev: CDR3β ni med celicami iz krvi ostalih bolnikov
        ucne = set(xl[(xl[:, c_p] != pid) & (xl[:, c_t] == 0), c_c].astype(int).tolist())
        nove = np.array([t not in ucne for t in tcr])

        zr = [r for r in zapisi if r['heldout_patient'] == pid]
        for r in zr:
            assert np.array_equal(np.asarray(r['x_true'], int), x.astype(int)), f'{r["_datoteka"]}: oznake se ne ujemajo'
        M = np.column_stack([pre] + [np.asarray(r['y_score'], float) for r in zr])
        imena = ['velikost_klona'] + [f'{r["variant"]}:s{r["seed"]}' + ('' if r['_glavni'] else ':ponovitev') for r in zr]
        auc, S = delong(x, M)
        B = boot_auc(x, M, tcr)
        auc_n = auc_stolpci(x[nove], M[nove])
        B_n = boot_auc(x[nove], M[nove], tcr[nove])
        n_info = {'n_celic_nove': int(nove.sum()), 'n_eksp_nove': int(x[nove].sum()),
                  'n_eksp_klonov_nove': int(len(np.unique(tcr[nove & x]))),
                  'n_celic': int(len(x)), 'n_eksp': int(x.sum()), 'n_eksp_klonov': int(len(np.unique(tcr[x])))}
        povz[f'bolnik_{pid}'] = n_info
        for i, ime in enumerate(imena):
            r = zr[i - 1] if i else None
            izh = max(r['baseline_results'].values()) if r and r.get('baseline_results') else np.nan
            se = np.sqrt(S[i, i])
            bl, bh = interval(B[:, i]); nl, nh = interval(B_n[:, i])
            v_zag.append({'bolnik': pid, 'ocena': ime, 'razlicica': r['variant'] if r else 'velikost_klona',
                          'seme': r['seed'] if r else None, 'glavni': bool(r['_glavni']) if r else True,
                          'auc': auc[i], 'delong_sp': max(0, auc[i] - Z95 * se), 'delong_zg': min(1, auc[i] + Z95 * se),
                          'boot_sp': bl, 'boot_zg': bh,
                          'auc_nove': auc_n[i], 'boot_nove_sp': nl, 'boot_nove_zg': nh,
                          'najboljse_izhodisce': izh,
                          'nad_izhodiscem_delong': bool(auc[i] - Z95 * se > izh) if r else None,
                          'nad_izhodiscem_boot': bool(bl > izh) if r else None,
                          'nad_izhodiscem_nove_boot': bool(nl > izh) if r else None})
        I = {ime: i for i, ime in enumerate(imena)}
        # različice do RNA + TCR pri istem semenu
        for seme in (0, 1, 2):
            a = I.get(f'rna_tcr:s{seme}')
            if a is None:
                continue
            for razl in ('rna_flux_tcr', 'flux_tcr'):
                b = I.get(f'{razl}:s{seme}')
                if b is None:
                    continue
                d = B[:, b] - B[:, a]; dl, dh = interval(d)
                var = S[a, a] + S[b, b] - 2 * S[a, b]
                v_razl.append({'bolnik': pid, 'seme': seme, 'razlicica': razl, 'delta_auc': auc[b] - auc[a],
                               'boot_sp': dl, 'boot_zg': dh, 'boot_p': p_boot(d),
                               'delong_p': 2 * norm.sf(abs((auc[b] - auc[a]) / np.sqrt(var)))})
        # semena iste različice (in ponovljeni zagon)
        for razl in RAZLICICE:
            k = [ime for ime in imena if ime.startswith(razl + ':') and not ime.endswith('ponovitev')]
            pari = list(combinations(k, 2))
            pon = [ime for ime in imena if ime.startswith(razl + ':') and ime.endswith('ponovitev')]
            pari += [(p.replace(':ponovitev', ''), p) for p in pon]
            for ia, ib in pari:
                a, b = I[ia], I[ib]
                d = B[:, b] - B[:, a]; dl, dh = interval(d)
                var = S[a, a] + S[b, b] - 2 * S[a, b]
                v_sem.append({'bolnik': pid, 'razlicica': razl, 'a': ia, 'b': ib, 'delta_auc': auc[b] - auc[a],
                              'boot_sp': dl, 'boot_zg': dh, 'boot_p': p_boot(d),
                              'delong_p': 2 * norm.sf(abs((auc[b] - auc[a]) / np.sqrt(var)))})

    zag = pd.DataFrame(v_zag); razl = pd.DataFrame(v_razl); sem = pd.DataFrame(v_sem)
    zag.to_csv(os.path.join(IZHOD, 'zagoni.csv'), index=False)
    razl.to_csv(os.path.join(IZHOD, 'razlicice.csv'), index=False)
    sem.to_csv(os.path.join(IZHOD, 'semena.csv'), index=False)
    glavni = sem[~sem.b.str.endswith('ponovitev')]
    povz['semena_znacilno_delong'] = int((glavni.delong_p < 0.05).sum())
    povz['semena_znacilno_boot'] = int((glavni.boot_p < 0.05).sum())
    povz['semena_pari'] = int(len(glavni))
    t = zag[zag.razlicica != 'velikost_klona']
    t = t[t.glavni]
    for razl_ in RAZLICICE:
        z = t[t.razlicica == razl_]
        povz[razl_] = {'nad_izhodiscem_delong': int(z.nad_izhodiscem_delong.sum()),
                       'nad_izhodiscem_boot': int(z.nad_izhodiscem_boot.sum()),
                       'nad_izhodiscem_nove_boot': int(z.nad_izhodiscem_nove_boot.sum()), 'zagonov': int(len(z))}
    json.dump(povz, open(os.path.join(IZHOD, 'povzetek.json'), 'w'), indent=2)
    pd.set_option('display.width', 250)
    print(zag.round(3).to_string())
    print(razl.round(3).to_string())
    print(sem.round(3).to_string())
    print(json.dumps(povz, indent=1))


if __name__ == '__main__':
    main()
