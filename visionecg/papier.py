"""Génération d'un tracé ECG de synthèse imprimé sur papier millimétré (image PNG).

Convention du papier : 25 mm/s et 10 mm/mV ; petits carreaux de 1 mm, gros de 5 mm ;
impulsion d'étalonnage de 1 mV pendant 0,2 s au début de la bande.
Le signal est fabriqué avec des ondes gaussiennes (P, Q, R, S, T) : aucune donnée réelle.
"""
import numpy as np
from PIL import Image, ImageDraw

VITESSE = 25.0  # mm/s
GAIN = 10.0     # mm/mV

# (position relative au pic R en s, amplitude en mV, largeur en s)
ONDES = {"P": (-0.20, 0.15, 0.025), "Q": (-0.03, -0.12, 0.010),
         "R": (0.0, 1.2, 0.012), "S": (0.03, -0.25, 0.012), "T": (0.30, 0.35, 0.045)}


def battements(duree, fc, irregularite=0.0, graine=0, debut=0.6):
    """Instants des pics R : rythme régulier, ou irrégulier (écart-type relatif des RR)."""
    rng = np.random.default_rng(graine)
    t, r = debut, []
    while t < duree - 0.4:
        r.append(t)
        rr = 60 / fc * (1 + irregularite * rng.standard_normal())
        t += max(0.3, rr)
    return np.array(r)


def signal(t, pics_r, largeur_qrs=1.0, ondes=ONDES):
    """Tension (mV) aux instants t. largeur_qrs > 1 élargit Q, R et S (bloc de branche)."""
    v = np.zeros_like(t, dtype=float)
    for tr in pics_r:
        for nom, (pos, amp, larg) in ondes.items():
            k = largeur_qrs if nom in "QRS" else 1.0
            v += amp * np.exp(-0.5 * ((t - tr - pos * (k if nom in "QS" else 1)) / (larg * k)) ** 2)
    return v


def imprimer(duree=10.0, fc=72, irregularite=0.0, largeur_qrs=1.0, px_par_mm=8, graine=0,
             hauteur_mm=40, bruit_mv=0.0):
    """Renvoie (image PIL, vérité terrain)."""
    rng = np.random.default_rng(graine + 1)
    marge_mm = 10  # place pour l'impulsion d'étalonnage
    largeur_mm = marge_mm + duree * VITESSE
    w, h = int(round(largeur_mm * px_par_mm)), int(round(hauteur_mm * px_par_mm))
    im = Image.new("RGB", (w, h), (255, 246, 246))
    d = ImageDraw.Draw(im)
    for i in range(int(largeur_mm) + 1):  # grille
        x = round(i * px_par_mm)
        d.line([(x, 0), (x, h)], fill=(240, 170, 170) if i % 5 == 0 else (250, 215, 215), width=2 if i % 5 == 0 else 1)
    for j in range(int(hauteur_mm) + 1):
        y = round(j * px_par_mm)
        d.line([(0, y), (w, y)], fill=(240, 170, 170) if j % 5 == 0 else (250, 215, 215), width=2 if j % 5 == 0 else 1)
    base = h * 0.6
    mv = px_par_mm * GAIN  # pixels par mV
    # impulsion d'étalonnage : 1 mV pendant 0,2 s, entre 2 et 7 mm
    x0, x1 = 2 * px_par_mm, (2 + 0.2 * VITESSE) * px_par_mm
    d.line([(0, base), (x0, base), (x0, base - mv), (x1, base - mv), (x1, base), (marge_mm * px_par_mm, base)],
           fill=(20, 20, 20), width=2)
    pics = battements(duree, fc, irregularite, graine)
    xs = np.arange(marge_mm * px_par_mm, w)
    t = (xs / px_par_mm - marge_mm) / VITESSE
    v = signal(t, pics, largeur_qrs) + bruit_mv * rng.standard_normal(len(t))
    ys = base - v * mv
    d.line(list(zip(xs.tolist(), ys.tolist())), fill=(20, 20, 20), width=2)
    verite = {"pics_r": pics, "fc": 60 / np.diff(pics).mean(), "px_par_mm": px_par_mm,
              "t": t, "v": v, "marge_mm": marge_mm, "amplitude_r": float(signal(np.array([pics[1]]), pics, largeur_qrs)[0]
                                     - np.median(signal(pics[1] - np.linspace(0.07, 0.12, 13), pics, largeur_qrs)))}
    return im, verite
