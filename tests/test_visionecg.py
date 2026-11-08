"""Vérifications : python tests/test_visionecg.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from visionecg import imprimer, analyser, numeriser, pas_grille, masques, charger

n_ok = 0


def ok(cond, msg):
    global n_ok
    if not cond:
        raise AssertionError(msg)
    n_ok += 1


# Grille : pixels par mm retrouvés, y compris pour des résolutions non entières
for px in (6, 8, 10.5, 12):
    im, _ = imprimer(duree=6, px_par_mm=px)
    est = pas_grille(masques(charger(im))[1])
    ok(abs(est - px) / px < 0.01, f"grille {px} px/mm -> {est:.3f}")

# Étalonnage 10 mm/mV retrouvé sur l'impulsion de 1 mV
for px in (6, 8, 12):
    _, _, info = numeriser(imprimer(duree=6, px_par_mm=px)[0])
    ok(abs(info["gain_mesure_mm_par_mv"] - 10) < 0.4, f"étalonnage ({px} px/mm) : {info['gain_mesure_mm_par_mv']:.2f}")

# Numérisation : fidèle au signal vrai en dehors des QRS (pentes trop raides pour un pixel)
for px in (8, 12):
    im, v = imprimer(duree=8, fc=70, px_par_mm=px, graine=4)
    t, vv, _ = numeriser(im)
    vt = np.interp(t, v["t"], v["v"])
    m = np.ones(len(t), bool)
    for r in v["pics_r"]:
        m &= np.abs(t - r) > 0.06
    rms = np.sqrt(np.mean((vt - vv)[m] ** 2))
    ok(rms < 0.02, f"erreur de numérisation hors QRS ({px} px/mm) : {rms:.4f} mV")

# Fréquence, battements, amplitude, sur plusieurs rythmes et tirages
for fc in (45, 60, 72, 100, 140):
    for g in range(2):
        im, v = imprimer(fc=fc, graine=g)
        r = analyser(im)
        ok(r["nb_battements"] == len(v["pics_r"]), f"FC {fc} : {r['nb_battements']} battements / {len(v['pics_r'])}")
        ok(abs(r["fc"] - v["fc"]) < 1.0, f"FC {fc} : {r['fc']:.1f}")
        ok(np.max(np.abs(r["instants_r"] - v["pics_r"])) < 0.012, f"FC {fc} : pics R à moins de 12 ms")
        ok(abs(r["amplitude_r_mv"] - v["amplitude_r"]) < 0.06, f"FC {fc} : amplitude R")
ok(analyser(imprimer(fc=45)[0])["frequence"] == "bradycardie", "bradycardie reconnue")
ok(analyser(imprimer(fc=130)[0])["frequence"] == "tachycardie", "tachycardie reconnue")

# Rythme irrégulier et variabilité
for g in range(3):
    im, v = imprimer(fc=90, irregularite=0.18, graine=10 + g)
    r = analyser(im)
    vrai_sdnn = np.diff(v["pics_r"]).std(ddof=1) * 1000
    ok(r["rythme"] == "irrégulier", f"rythme irrégulier reconnu (tirage {g})")
    ok(abs(r["sdnn_ms"] - vrai_sdnn) < 8, f"SDNN {r['sdnn_ms']:.0f} / {vrai_sdnn:.0f} ms")
ok(analyser(imprimer(fc=75, graine=3)[0])["rythme"] == "régulier", "rythme régulier reconnu")

# QRS : élargi détecté, et durée croissante avec la largeur des ondes
q = [analyser(imprimer(fc=65, largeur_qrs=k)[0])["qrs_ms"] for k in (1.0, 1.6, 2.4)]
ok(q[0] < q[1] < q[2], f"durée du QRS croissante {q}")
ok(q[0] < 110 and q[2] >= 120, "QRS fin / QRS élargi")

# Robustesse : un peu de bruit sur le tracé
im, v = imprimer(fc=80, bruit_mv=0.02, graine=6)
r = analyser(im)
ok(r["nb_battements"] == len(v["pics_r"]) and abs(r["fc"] - v["fc"]) < 1.5, "tracé bruité")

print(f"{n_ok} vérifications réussies")
