# -*- coding: utf-8 -*-
"""
Enotno belezenje eksperimentov TRIM-Flux -> en JSON na zagon (na Drive).

Uporaba (v eval notebooku, na koncu):
    from experiment_logger import log_run   # ali prilepi funkcijo v celico
    log_run(
        log_dir=os.path.join(data_parent_folder, 'results/logs'),
        variant='flux_tcr',            # 'rna_tcr' | 'flux_tcr' | 'rna_flux_tcr'
        heldout_patient=heldout_patient,
        seed=seed,                     # ce ga imas; sicer None
        roc_auc=roc_auc,
        x=x, y=y,                      # resnica (0/1) in napoved (zvezna)
        thresholds=threshold_reports,  # dict {ime: pred_array} ali None -> sam izracuna
        baseline_results=baseline_results if 'baseline_results' in dir() else None,
        config=config if 'config' in dir() else None,     # hiperparametri iz treninga
        thresh_fitted=thresh_fitted if 'thresh_fitted' in dir() else None,
        extra={'n_train_cells': int(mask_train.sum()) if 'mask_train' in dir() else None},
        timestamp=RUN_TS,              # PODAJ iz notebooka! (time.strftime), glej opombo spodaj
    )

OPOMBA o timestampu: podaj ga iz notebooka (npr. RUN_TS = time.strftime('%Y%m%d-%H%M%S')),
NE znotraj te funkcije z Date.now ipd. Tako je ime datoteke stabilno in znano ze pred klicem.
"""
import os
import json
import numpy as np


def _confusion_at(x, y, thresh):
    """Vrne confusion + metrike pri danem pragu (brez sklearn odvisnosti — cist numpy)."""
    pred = (y >= thresh).astype(int) if thresh > 0 else (y > 0).astype(int)
    x = np.asarray(x).astype(int)
    tp = int(((pred == 1) & (x == 1)).sum())
    tn = int(((pred == 0) & (x == 0)).sum())
    fp = int(((pred == 1) & (x == 0)).sum())
    fn = int(((pred == 0) & (x == 1)).sum())
    def safe(n, d):
        return float(n / d) if d else 0.0
    prec_exp = safe(tp, tp + fp)
    rec_exp  = safe(tp, tp + fn)
    f1_exp   = safe(2 * prec_exp * rec_exp, prec_exp + rec_exp)
    prec_not = safe(tn, tn + fn)
    rec_not  = safe(tn, tn + fp)
    f1_not   = safe(2 * prec_not * rec_not, prec_not + rec_not)
    return {
        'thresh': float(thresh),
        'n_pred_expanded': int(pred.sum()),
        'confusion': {'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp},
        'EXPANDED': {'precision': round(prec_exp, 4), 'recall': round(rec_exp, 4), 'f1': round(f1_exp, 4)},
        'NOT':      {'precision': round(prec_not, 4), 'recall': round(rec_not, 4), 'f1': round(f1_not, 4)},
        'specificity': round(safe(tn, tn + fp), 4),
    }


def _to_jsonable(o):
    """Rekurzivno pretvori numpy tipe v cist python (JSON ne pozna np.int64/float32)."""
    if isinstance(o, dict):
        return {str(k): _to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_to_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def log_run(log_dir, variant, heldout_patient, roc_auc, x, y,
            timestamp, seed=None, baseline_results=None, config=None,
            thresh_fitted=None, extra=None):
    """
    Zabelezi en zagon kot JSON. Vrne pot do shranjene datoteke.
    Robustno: napaka pri enem polju ne sme ubiti eksperimenta -> sglobal try/except.
    """
    try:
        os.makedirs(log_dir, exist_ok=True)
        x = np.asarray(x)
        y = np.asarray(y)
        pid = heldout_patient[0] if isinstance(heldout_patient, (list, tuple)) else int(heldout_patient)

        # --- pragovi: y>0, Youden, enaka bazna stopnja ---
        thresholds = {}
        n_pos = int((x == 1).sum())
        both_classes = 0 < n_pos < len(x)
        if both_classes:
            thresholds['y>0'] = _confusion_at(x, y, 0.0)
            # Youden (max TPR-FPR) — sam izracuna prag brez sklearn
            order = np.argsort(-y)
            xs = x[order]
            P, N = n_pos, len(x) - n_pos
            tp_c = np.cumsum(xs == 1)
            fp_c = np.cumsum(xs == 0)
            tpr = tp_c / P if P else np.zeros_like(tp_c, dtype=float)
            fpr = fp_c / N if N else np.zeros_like(fp_c, dtype=float)
            j_idx = int(np.argmax(tpr - fpr))
            youden_thresh = float(y[order][j_idx])
            thresholds['youden'] = _confusion_at(x, y, max(youden_thresh, 1e-12))
            base_rate = float(x.mean())
            q_thresh = float(np.quantile(y, 1 - base_rate))
            thresholds['base_rate'] = _confusion_at(x, y, max(q_thresh, 1e-9))

        record = {
            'variant': variant,                       # rna_tcr | flux_tcr | rna_flux_tcr
            'heldout_patient': pid,
            'seed': seed,
            'timestamp': timestamp,
            'roc_auc': round(float(roc_auc), 4) if roc_auc is not None else None,
            'n_cells_eval': int(len(x)),
            'n_expanded_true': n_pos,
            'expanded_frac': round(float(x.mean()), 4) if len(x) else None,
            'both_classes': bool(both_classes),
            'thresholds': thresholds,
            'baseline_results': _to_jsonable(baseline_results) if baseline_results else None,
            'thresh_fitted': float(thresh_fitted) if thresh_fitted is not None else None,
            'config': _to_jsonable(config) if config else None,
            'extra': _to_jsonable(extra) if extra else None,
            # surova x (resnica 0/1) + y (napoved) za POOLED ROC cez vec zagonov/pacientov
            'x_true': [int(v) for v in x.tolist()],
            'y_score': [float(v) for v in y.tolist()],
        }

        seed_str = f'_seed{seed}' if seed is not None else ''
        fname = f'{variant}_P{pid}{seed_str}_{timestamp}.json'
        path = os.path.join(log_dir, fname)
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(_to_jsonable(record), f, ensure_ascii=False, indent=2)
        print(f'[log_run] shranjeno -> {path}')
        print(f'          variant={variant} P{pid} seed={seed} AUC={record["roc_auc"]}')
        return path
    except Exception as e:
        # log NE sme ubiti eksperimenta
        print(f'[log_run] OPOZORILO: belezenje ni uspelo ({e}). Rezultat NI izgubljen, le ni logiran.')
        return None


def collect_logs(log_dir):
    """Zbere vse JSON loge v en pandas DataFrame (za primerjavo cez zagone/paciente/variante)."""
    import pandas as pd
    rows = []
    for fn in sorted(os.listdir(log_dir)):
        if not fn.endswith('.json'):
            continue
        try:
            r = json.load(open(os.path.join(log_dir, fn), encoding='utf-8'))
        except Exception:
            continue
        row = {
            'variant': r.get('variant'), 'patient': r.get('heldout_patient'),
            'seed': r.get('seed'), 'auc': r.get('roc_auc'),
            'n_expanded': r.get('n_expanded_true'), 'timestamp': r.get('timestamp'),
        }
        yo = (r.get('thresholds') or {}).get('youden') or {}
        row['exp_recall_youden'] = yo.get('EXPANDED', {}).get('recall')
        row['exp_precision_youden'] = yo.get('EXPANDED', {}).get('precision')
        row['exp_f1_youden'] = yo.get('EXPANDED', {}).get('f1')
        rows.append(row)
    df = pd.DataFrame(rows)
    return df


def summarize_by_variant(log_dir):
    """Povzetek: mean +/- std AUC po (variant, patient) cez seede — za Monte Carlo primerjavo."""
    import pandas as pd
    df = collect_logs(log_dir)
    if df.empty:
        print('Ni logov.')
        return df
    g = df.groupby(['variant', 'patient'])['auc'].agg(['mean', 'std', 'count']).reset_index()
    g['mean'] = g['mean'].round(4)
    g['std'] = g['std'].round(4)
    print(g.to_string(index=False))
    return g


def pooled_roc(log_dir, variant=None, seed=None):
    """
    POOLED ROC: zdruzi surove x/y cez VSE zagone (paciente/seede) v en ROC -> ozji CI.
    Kot original (vse heldout celice v en pool). Vrne dict {variant: (auc, n_cells, n_pos)}.
    Filtriraj z variant/seed za ozji nabor. Zahteva log-e z 'x_true'/'y_score'
    (novi log_run jih shrani; stari log-i brez tega se preskocijo z opozorilom).
    """
    import numpy as np
    from collections import defaultdict
    pools = defaultdict(lambda: {'x': [], 'y': [], 'n_runs': 0})
    skipped = 0
    for fn in sorted(os.listdir(log_dir)):
        if not fn.endswith('.json'):
            continue
        try:
            r = json.load(open(os.path.join(log_dir, fn), encoding='utf-8'))
        except Exception:
            continue
        if variant is not None and r.get('variant') != variant:
            continue
        if seed is not None and r.get('seed') != seed:
            continue
        xt, ys = r.get('x_true'), r.get('y_score')
        if xt is None or ys is None:
            skipped += 1
            continue
        v = r.get('variant')
        pools[v]['x'].extend(xt)
        pools[v]['y'].extend(ys)
        pools[v]['n_runs'] += 1
    if skipped:
        print(f'[pooled_roc] {skipped} starih logov brez x_true/y_score preskocenih.')

    out = {}
    try:
        import sklearn.metrics
        have_skl = True
    except Exception:
        have_skl = False
    for v, d in pools.items():
        x = np.asarray(d['x']); y = np.asarray(d['y'])
        n_pos = int((x == 1).sum())
        if not (0 < n_pos < len(x)):
            print(f'[pooled_roc] {v}: en sam razred ({n_pos}/{len(x)}) -> ROC ni definiran.')
            out[v] = (None, len(x), n_pos)
            continue
        if have_skl:
            fpr, tpr, _ = sklearn.metrics.roc_curve(x, y, pos_label=1)
            auc = float(sklearn.metrics.auc(fpr, tpr))
        else:  # cist numpy AUC (rank)
            order = np.argsort(-y); xs = x[order]
            P, N = n_pos, len(x) - n_pos
            tpr = np.cumsum(xs == 1) / P; fpr = np.cumsum(xs == 0) / N
            auc = float(np.trapz(np.r_[0, tpr], np.r_[0, fpr]))
        out[v] = (round(auc, 4), len(x), n_pos)
        print(f'[pooled_roc] {v}: AUC={auc:.4f}  (pooled {d["n_runs"]} zagonov, {len(x)} celic, {n_pos} ekspandiranih)')
    return out
