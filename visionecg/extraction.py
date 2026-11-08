"""Extraction de paramètres physiologiques à partir d'une image d'ECG sur papier.

Étapes :
1. séparer le tracé (pixels sombres) de la grille (pixels rosés) ;
2. mesurer le pas de la grille (autocorrélation du profil de la grille) : pixels par mm ;
3. mesurer l'impulsion d'étalonnage (1 mV) pour vérifier le gain ;
4. numériser : une valeur de tension par colonne de pixels, convertie en mV et en secondes ;
5. détecter les pics R, puis calculer FC, intervalles RR, variabilité, amplitude R, durée du QRS.
"""
import numpy as np
from PIL import Image

VITESSE = 25.0
GAIN = 10.0


def charger(image):
    if isinstance(image, Image.Image):
        return np.asarray(image.convert("RGB"), float)
    return np.asarray(Image.open(image).convert("RGB"), float)


def masques(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    trace = (r < 110) & (g < 110) & (b < 110)
    grille = (r > 200) & (g < 235) & (r - g > 15) & ~trace
    return trace, grille


def pas_grille(grille):
    """Pixels par mm : période fondamentale du profil des lignes verticales de la grille."""
    profil = grille.sum(0).astype(float)
    profil -= profil.mean()
    n = len(profil)
    ac = np.correlate(profil, profil, mode="full")[n - 1:]
    ac /= ac[0]
    debut = 3
    candidats = [k for k in range(debut, min(200, n // 4)) if ac[k] > ac[k - 1] and ac[k] >= ac[k + 1] and ac[k] > 0.1]
    if not candidats:
        raise ValueError("grille introuvable")
    k = candidats[0]
    # affiner avec les gros carreaux (5 mm), plus précis
    k5 = int(round(5 * k))
    fen = range(max(k5 - 3, 1), k5 + 4)
    k5 = max(fen, key=lambda i: ac[i])
    a, b, c = ac[k5 - 1], ac[k5], ac[k5 + 1]
    dec = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) else 0
    return (k5 + dec) / 5


def courbe(trace, ecart_min=3):
    """Ligne du tracé pour chaque colonne, NaN si absente.

    Sur les parties plates : moyenne des pixels sombres de la colonne. Sur les ondes
    raides (pic R, onde S) la ligne est presque verticale et une colonne contient tout le
    segment : on garde alors le pixel le plus éloigné de la ligne de base, sinon les pics
    seraient sous-estimés.
    """
    h, w = trace.shape
    lignes = np.arange(h)[:, None]
    n = trace.sum(0)
    with np.errstate(invalid="ignore"):
        moy = (trace * lignes).sum(0) / n
    moy[n == 0] = np.nan
    base = np.nanmedian(moy)
    haut = np.where(trace, lignes, h).min(0).astype(float)
    bas = np.where(trace, lignes, -1).max(0).astype(float)
    y = moy.copy()
    au_dessus = moy < base - ecart_min
    en_dessous = moy > base + ecart_min
    y[au_dessus] = haut[au_dessus] + 0.5  # le trait fait 2 px : on vise son milieu
    y[en_dessous] = bas[en_dessous] - 0.5
    y[n == 0] = np.nan
    return y


def etalonnage(y, px_mm):
    """Hauteur (px) de l'impulsion d'étalonnage dans les 10 premiers mm."""
    zone = y[: int(10 * px_mm)]
    zone = zone[~np.isnan(zone)]
    base = np.percentile(zone, 90)  # la ligne de base est plus bas sur l'image
    haut = np.percentile(zone, 10)
    return base - haut


def numeriser(image):
    """Renvoie (t en s, tension en mV, infos d'étalonnage)."""
    rgb = charger(image)
    trace, grille = masques(rgb)
    px_mm = pas_grille(grille)
    y = courbe(trace)
    px_mv_mesure = etalonnage(y, px_mm)
    debut = int(round(10 * px_mm))
    y = y[debut:]
    ok = ~np.isnan(y)
    y = np.interp(np.arange(len(y)), np.flatnonzero(ok), y[ok])
    base = np.median(y)
    px_mv = px_mm * GAIN
    v = (base - y) / px_mv
    t = np.arange(len(y)) / px_mm / VITESSE
    return t, v, {"px_par_mm": px_mm, "gain_mesure_mm_par_mv": px_mv_mesure / px_mm}


def reechantillonner(t, v, fe=250.0):
    tt = np.arange(0, t[-1], 1 / fe)
    return tt, np.interp(tt, t, v)


def detecter_r(t, v, fe=250.0):
    """Pics R : énergie de la dérivée, seuil adaptatif, période réfractaire de 250 ms."""
    d = np.gradient(v) * fe
    e = d ** 2
    k = int(0.1 * fe)
    e = np.convolve(e, np.ones(k) / k, mode="same")
    seuil = 0.3 * np.percentile(e, 99)
    pics = []
    i = 0
    refr = int(0.25 * fe)
    while i < len(e):
        if e[i] > seuil:
            fin = min(len(v), i + int(0.15 * fe))
            debut = max(0, i - int(0.05 * fe))
            j = debut + int(np.argmax(v[debut:fin]))
            if not pics or j - pics[-1] > refr:
                pics.append(j)
            i = j + refr
        else:
            i += 1
    return np.array(pics)


def duree_qrs(t, v, pics, fe=250.0):
    """Durée du QRS : autour de chaque R, zone où la pente dépasse 15 % de la pente maximale."""
    d = np.abs(np.gradient(v) * fe)
    durees = []
    for p in pics:
        a, b = max(0, p - int(0.12 * fe)), min(len(v), p + int(0.12 * fe))
        seg = d[a:b]
        if len(seg) < 5:
            continue
        idx = np.flatnonzero(seg > 0.15 * seg.max())
        durees.append((idx[-1] - idx[0]) / fe)
    return float(np.median(durees)) if durees else float("nan")


def amplitude_r(t, v, pics, fe=250.0):
    """Hauteur du pic R au-dessus de la ligne isoélectrique, mesurée sur le segment PR
    (de 120 à 70 ms avant R) : plus juste que la médiane du signal quand le cœur bat vite."""
    h = []
    for p in pics:
        a, b = p - int(0.12 * fe), p - int(0.07 * fe)
        if a >= 0:
            h.append(v[p] - np.median(v[a:b]))
    return float(np.median(h)) if h else float("nan")


def analyser(image):
    t, v, info = numeriser(image)
    tt, vv = reechantillonner(t, v)
    pics = detecter_r(tt, vv)
    rr = np.diff(tt[pics])
    res = dict(info)
    res.update({
        "nb_battements": len(pics),
        "instants_r": tt[pics],
        "fc": 60 / rr.mean() if len(rr) else float("nan"),
        "rr_moyen_ms": rr.mean() * 1000 if len(rr) else float("nan"),
        "sdnn_ms": rr.std(ddof=1) * 1000 if len(rr) > 2 else float("nan"),
        "amplitude_r_mv": amplitude_r(tt, vv, pics),
        "qrs_ms": duree_qrs(tt, vv, pics) * 1000,
        "signal": (tt, vv),
    })
    cv = rr.std() / rr.mean() if len(rr) > 2 else 0
    res["rythme"] = "irrégulier" if cv > 0.08 else "régulier"
    fc = res["fc"]
    res["frequence"] = "bradycardie" if fc < 60 else "tachycardie" if fc > 100 else "normale"
    res["qrs"] = "élargi (≥ 120 ms)" if res["qrs_ms"] >= 120 else "fin"
    return res
