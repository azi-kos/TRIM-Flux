"""Od kod pride AUC modela TRIM: delitev na y = 0 in y > 0 ter združevanje semen.

Bere zapise zagonov v results/logs/logs (x_true, y_score) in za vsakega bolnika
in različico izračuna:
- AUC posameznega zagona,
- AUC samo iz delitve y > 0 (celica je pred zdravljenjem v psevdoklonu z vsaj
  dvema celicama),
- AUC znotraj celic z y > 0,
- AUC povprečja ocen in povprečja rangov po semenih,
- AUC deleža zagonov, v katerih ima celica y > 0,
- Jaccardov indeks množic celic z y > 0 med zagoni.

Pri podvojenem zagonu (flux_tcr, bolnik 23, seme 0) uporabi prvega, kot
rezultati_napovedi.py. Izhod: diploma/rezultati/psevdokloni_auc.csv in
diploma/rezultati/psevdokloni_jaccard.csv.

Zagon iz korena repozitorija: python3 diploma/skripte/analiza_psevdoklonov.py
"""
import csv
import glob
import json
import os
from collections import defaultdict

import numpy as np

KOREN = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ZAPISI = os.path.join(KOREN, 'results', 'logs', 'logs')
IZHOD = os.path.join(KOREN, 'diploma', 'rezultati')
RAZLICICE = ('rna_tcr', 'flux_tcr', 'rna_flux_tcr')


def auc(x, s):
    """AUC kot verjetnost, da ima pozitivna celica višjo oceno od negativne (vezi 1/2)."""
    x = np.asarray(x)
    s = np.asarray(s, float)
    p, n = s[x == 1], s[x == 0]
    if len(p) == 0 or len(n) == 0:
        return float('nan')
    return float(((p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum()) / (len(p) * len(n)))


def rangi_s_vezmi(y):
    """Povprečni rangi (vezi dobijo povprečje svojih rangov), kot scipy.stats.rankdata."""
    y = np.asarray(y, float)
    red = np.argsort(y, kind='mergesort')
    ys = y[red]
    r = np.empty(len(y))
    i = 0
    while i < len(y):
        j = i
        while j + 1 < len(y) and ys[j + 1] == ys[i]:
            j += 1
        r[red[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return r


def nalozi():
    zagoni = defaultdict(list)
    for pot in sorted(glob.glob(os.path.join(ZAPISI, '*.json'))):
        with open(pot, encoding='utf-8') as f:
            r = json.load(f)
        zagoni[(r['variant'], int(r['heldout_patient']), int(r['seed']))].append(
            (r['timestamp'], np.array(r['x_true']), np.array(r['y_score'], float)))
    return {k: sorted(v, key=lambda t: t[0])[0] for k, v in zagoni.items()}


def main():
    glavni = nalozi()
    vrstice, jaccard = [], []
    for bolnik in sorted({k[1] for k in glavni}):
        kljuci = sorted(k for k in glavni if k[1] == bolnik)
        x = glavni[kljuci[0]][1]
        assert all(np.array_equal(x, glavni[k][1]) for k in kljuci), 'oznake se med zagoni razlikujejo'
        poz = {k: glavni[k][2] > 0 for k in kljuci}
        for razl in RAZLICICE + ('vse',):
            kr = [k for k in kljuci if razl == 'vse' or k[0] == razl]
            Y = np.array([glavni[k][2] for k in kr])
            rangi = np.array([rangi_s_vezmi(y) for y in Y])
            for k, y in zip(kr, Y):
                if razl == 'vse':
                    break
                m = y > 0
                vrstice.append({
                    'bolnik': bolnik, 'razlicica': razl, 'seme': k[2],
                    'n_celic': len(x), 'n_ekspandiranih': int(x.sum()),
                    'auc': round(auc(x, y), 4),
                    'auc_samo_y_poz': round(auc(x, m.astype(float)), 4),
                    'auc_znotraj_y_poz': round(auc(x[m], y[m]), 4),
                    'delez_y0': round(float((~m).mean()), 4),
                    'delez_eksp_y0': round(float(x[~m].mean()), 4) if (~m).any() else '',
                    'delez_eksp_y_poz': round(float(x[m].mean()), 4) if m.any() else '',
                })
            vrstice.append({
                'bolnik': bolnik, 'razlicica': razl, 'seme': 'zdruzeno',
                'n_celic': len(x), 'n_ekspandiranih': int(x.sum()),
                'auc': round(auc(x, Y.mean(0)), 4),
                'auc_povprecje_rangov': round(auc(x, rangi.mean(0)), 4),
                'auc_delez_zagonov_y_poz': round(auc(x, np.mean([poz[k] for k in kr], 0)), 4),
                'n_zagonov': len(kr),
            })
        for a in kljuci:
            for b in kljuci:
                if a < b:
                    jaccard.append({
                        'bolnik': bolnik, 'zagon_a': f'{a[0]}:{a[2]}', 'zagon_b': f'{b[0]}:{b[2]}',
                        'jaccard_y_poz': round(float((poz[a] & poz[b]).sum() / max(1, (poz[a] | poz[b]).sum())), 4),
                    })
    os.makedirs(IZHOD, exist_ok=True)
    polja = ['bolnik', 'razlicica', 'seme', 'n_celic', 'n_ekspandiranih', 'auc', 'auc_samo_y_poz',
             'auc_znotraj_y_poz', 'delez_y0', 'delez_eksp_y0', 'delez_eksp_y_poz',
             'auc_povprecje_rangov', 'auc_delez_zagonov_y_poz', 'n_zagonov']
    with open(os.path.join(IZHOD, 'psevdokloni_auc.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=polja)
        w.writeheader()
        w.writerows(vrstice)
    with open(os.path.join(IZHOD, 'psevdokloni_jaccard.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['bolnik', 'zagon_a', 'zagon_b', 'jaccard_y_poz'])
        w.writeheader()
        w.writerows(jaccard)
    for v in vrstice:
        print({k: v[k] for k in v if v[k] != ''})


if __name__ == '__main__':
    main()
